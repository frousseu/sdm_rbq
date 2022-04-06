

library(atlasBE)
library(dbplyr)
library(dplyr)
library(magrittr)
library(DBI)
library(sf)
library(data.table)
library(foreach)
library(doParallel)
library(raster)
library(mapSpecies)
library(concaveman)
library(scales)
library(RPostgreSQL)
library(rpostgis)

#load("/data/predictors_sdm/predictors.RData")

login <- readRDS("/data/atlasDB/.auth")
con <- conn(login$user, login$pwd, login$host, login$dbname)

# works, but says non-finite coordinates
pgListGeom(con, geog = TRUE)
s<-pgGetGeom(
  con,
  name=c("public","observations"),
  geom = "geom",
  other.cols = c("id_taxa"),
  query = NULL
)

# runs for 30 min + and then uses all memory and crashes
# s<-st_read(con, "observations", query = "select id,geom from observations limit 3;")



DBI::dbListTables(con)
dbListFields(con,"observations")

sp<-dplyr::tbl(con, "taxa") %>% collect()
taxonomy<-dplyr::tbl(con,"public.taxa") %>% collect()
variables<-dplyr::tbl(con,"variables") %>% collect()

d<-dplyr::tbl(con,"observations") %>% 
   dplyr::left_join(dplyr::tbl(con,"taxa"),by=c("id_taxa"="id")) %>% 
   dplyr::left_join(dplyr::tbl(con,"variables"),by=c("id_variables"="id")) %>% 
   #dplyr::select(id_taxa,scientific_name,rank,species_gr,id_variables,name,obs_value,geom) %>% 
   dplyr::select(scientific_name,rank,species_gr,geom) %>% 
   dplyr::filter(rank=="species") %>%
   #dplyr::filter(species_gr%in%c("Angiospermes","Conifères","Bryophytes")) %>%
   #dplyr::filter(species_gr%in%c("Amphibiens","Reptiles")) %>%
   dplyr::filter(species_gr%in%c("Poissons")) %>%
   dplyr::select(scientific_name,geom) %>% 
   dplyr::collect()

d[,c("lon", "lat")]<-sf::st_coordinates(sf::st_as_sfc(structure(d$geom,class='WKB'),EWKB=TRUE,crs=4326))
d<-d %>% dplyr::select(-geom)
d$species<-d$scientific_name

s<-st_as_sf(d,coords=c("lon","lat"),crs=4326)
s<-st_transform(s,prj)

e<-rasterize(s,predictors[[1]],field=1,fun="count",background=0)
names(e)<-"sbias"
rp<-stack(predictors,e)

region<-st_buffer(concaveman(s,concavity=2),100)
region<-as(region,"Spatial")
region<-spTransform(region,crs(rp)) ### added after for crs problem not sure if it works
#domain<-st_coordinates(st_cast(region,"MULTIPOINT"))[,1:2]
domain<-st_coordinates(st_cast(st_as_sf(region),"MULTIPOINT"))[,1:2]

pedge<-0.03
edge<-min(c(diff(bbox(region)[0.2,])*pedge,diff(bbox(region)[2,])*pedge))
edge


vars<-c("tmax","tmax2","forested","builtup","prec","prec2")
#vars<-names(predictors)
r<-stack(rp[[names(rp)[names(rp)%in%c(vars,"sbias")]]])

bpriors<-list(prec=list(default=1/(10)^2,latitude=1/(10)^2,latitude2=1/(10)^2,Intercept=1/(100)^2,sbias=1/(100)^2),mean=list(default=0,Intercept=0,latitude=0,latitude2=0,sbias=0))

Mesh <- inla.mesh.2d(loc.domain = domain, #coordinates(occsp),
                     max.edge = c(edge,edge*1.5),
                     min.angle = 21,
                     cutoff = edge/2,
                     offset = c(edge,edge*2),
                     crs = crs(region))

explana<-explanaMesh(sPoly=region,mesh=Mesh,X=r)

weight <- ppWeight(sPoly = region, mesh = Mesh)

f<-formula(paste("y~",paste(vars,collapse="+")))

#species<-c("Acer saccharum","Acer rubrum")
species<-names(rev(sort(table(s$species))))[1:10]
#f<-function(i){as.integer(any(i==sp))}
init<-vector(mode="list",length=length(species)*2)
ninit<-1

loccs<-lapply(species,function(x){
  occ<-s[s$species==x,]  
  occsp<-as(occ,"Spatial")
})
names(loccs)<-species

cl<-makeCluster(10)
registerDoParallel(cl)
foreach(j=1:length(loccs),
  .packages=c("INLA","mapSpecies","raster","sp","exactextractr","scales")
  #.export=c("loccs","explana","weight","bpriors","f","r","q","predictors")
) %dopar% {
########################
########################
sp<-names(loccs)[j]
#occ<-s[s$species==sp,]  
occsp<-loccs[[j]]
occ<-st_as_sf(occsp)
#explana<-explana
#weight<-weight
#bpriors<-bpriors
  
mpp <- ppSpace(f, sPoints = occsp,
               explanaMesh = explana,
               ppWeight = weight,
               prior.range = c(100,0.01),
               prior.sigma = c(0.1,0.01),
               num.threads=1,
               many=TRUE,
               control.inla = list(int.strategy = "eb"),
               fix = NULL,
               sboffset = "sbias",
               control.fixed = bpriors
)

colo<-colorRampPalette(c("grey90","steelblue","steelblue2","gold2","tomato2","red4"))(200)

png(paste0("/data/sdm_rbq/plots/",gsub(" ","_",sp),"_atlas.png"),units="in",width=10,height=10,res=300)
m<-list(mpp=mpp)
par(mfrow=n2mfrow(length(m)),mar=c(3,3,3,7))
for(i in 1:length(m)){
  mapMean<-mapSpace(m[[i]],dims=1*dim(r)[1:2],type="mean",sPoly=as(q,"Spatial"))
  #mapMean <- mask(mapMean,q[q$NAME_1=="Québec",]) #spacePoly
  mapMean <- mask(mapMean,q) #spacePoly
  buff<-st_as_sf(st_buffer(st_convex_hull(st_union(occ)),dist=500))
  mapMean <- mask(mapMean,buff) #spacePoly
  water<-predictors[["water"]]
  water[water<0.99]<-NA
  water<-resample(water,mapMean)
  mapMean <- mask(mapMean,water,inverse=TRUE) #spacePoly
  names(mapMean)<-paste(sp,names(m)[i])
  init[[ninit]]<-mapMean
  ninit<-ninit+1
  plot(buff,border=FALSE)
  plot(mapMean,add=TRUE, col = colo, axes = TRUE, box = FALSE, main = paste(names(m)[i],sp,sep=" - "))#,xlim=c(-1100,-900),ylim=c(200,300))
  #points(occsp,pch=16,col=alpha("black",0.95),cex=0.6)
  plot(st_geometry(q),add=TRUE,border=gray(0,0.15))
  ### contour
  ID<-inla.stack.index(attributes(m[[i]])$Stack,tag="est")$data
  #contour(mapMean,levels=quantile(m[[i]]$summary.fitted.values[["mean"]][ID],0.05),add=TRUE,col="red")
  lim<-rasterToContour(mapMean,levels=max(values(mapMean),na.rm=TRUE)*0.07)
  plot(lim,add=TRUE,col=alpha("steelblue",0.35),lwd=2)
  crs(mapMean)<-crs(occsp)
  writeRaster(mapMean, filename=paste0("/data/sdm_rbq/rasters/",gsub(" ","_",sp),"_atlas.tif"), format="GTiff", overwrite=TRUE)
}
par(mfrow=c(1,1))
dev.off()

print(j)
}






png("/data/sdm_rbq/plots/crap.png",units="in",res=200,width=7,height=7)
plot(Mesh)
plot(st_geometry(q),add=TRUE)
points(loccs[[6]])
dev.off()
