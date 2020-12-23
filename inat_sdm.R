

library(mapSpecies)
library(sf)
library(data.table)
library(scales)
library(raster)
library(concaveman)
library(fields)

#d<-fread("C:/Users/rouf1703/Documents/UdeS/Programmation/mapSpecies/observations-96473_20200626QC.csv")
#spinat<-d[,.N,by=.(taxon_species_name)][order(-N),][N>100,]$taxon_species_name
#d<-d[d$taxon_species_name%in%spinat,]
#d[,.(n=.N),by=.(taxon_species_name)][order(-n)][1:100]

keep<-c("kingdom","phylum","class","order","family","genus","species","countryCode","stateProvince","decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters","day","month","year","recordedBy")
d<-fread("D:/iNatGBIF/0137377-200613084148143.csv",encoding="UTF-8",select=keep)
d<-d[!is.na(d$decimalLatitude),]
d<-d[d$countryCode%in%c("US","CA"),]
prov<-c("Québec","Ontario","New Brunswick","Newfoundland and Labrador","Vermont","Maine","New Hampshire","New York","Massachusetts","Nunavut","Nova Scotia")
#d<-d[d$stateProvince%in%prov,]
d<-d[which(d$coordinateUncertaintyInMeters<50000),]
s<-st_as_sf(d,coords=c("decimalLongitude","decimalLatitude"),crs=4326)
s<-st_transform(s,prj)
#o<-st_intersects(s,q)

g<-st_make_grid(q,cellsize=10)
o<-st_intersects(s,g)
s$cell<-unlist(lapply(o,"[",1))
s<-s[!is.na(s$cell),]
s<-s[!duplicated(as.data.frame(s[,c("recordedBy","species","cell")])[,1:3]),]

#plot(st_geometry(q),border="grey85")
#plot(st_geometry(s),pch=16,col=gray(0,0.1),cex=0.5,add=TRUE)
#plot(log(r),axes=TRUE)
#plot(log(r),xlim=c(-400,800),ylim=c(-250,250))
#plot(st_geometry(q),border="grey65",add=TRUE)

#e<-rasterize(s,predictors[[1]],field=1,fun="count",background=0)[[1]]
e<-rasterize(s[s$class%in%c("Aves"),1],predictors[[1]],field=1,fun="count",background=0)
names(e)<-"sbias"
rp<-stack(predictors,e)#/10000)

#e2<-rasterize(occ,predictors[[1]],field=1,fun="count",background=0)

#################################################
#################################################


species<-c("Setophaga americana")
f<-function(i){as.integer(any(i==sp))}

init<-vector(mode="list",length=length(species)*2)
ninit<-1

for(j in 1:length(species)){
########################
########################
sp<-species[j]
occ<-s[s$species==sp,]  
occsp<-as(occ,"Spatial")
  
region<-st_buffer(concaveman(st_cast(q,"MULTIPOINT"),concavity=2),100)
domain<-st_coordinates(st_cast(region,"MULTIPOINT"))[,1:2]
region<-as(region,"Spatial")

pedge<-0.04
edge<-min(c(diff(bbox(region)[0.2,])*pedge,diff(bbox(region)[2,])*pedge))
edge

#vars<-c("tmean","tmean2","latitude","latitude2","prec","builtup","cultivated","conifers","forested")
#vars<-c("tmax","builtup")
#vars<-c("1")
#vars<-c("tmax","tmax2","trange","builtup","forested")
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
mpp <- ppSpace(f, sPoints = occsp,
               explanaMesh = explana,
               ppWeight = weight,
               prior.range = c(50,0.1),
               prior.sigma = c(0.5,0.1),
               num.threads=7,
               many=TRUE,
               control.inla = list(int.strategy = "eb"),
               fix = NULL,
               sboffset = "sbias",
               control.fixed = bpriors
)

colo<-colorRampPalette(c("grey90","steelblue","steelblue2","gold2","tomato2","red4"))(200)

m<-list(mpp=mpp)
par(mfrow=n2mfrow(length(m)),mar=c(3,3,3,7))
for(i in 1:length(m)){
  mapMean<-mapSpace(m[[i]],dims=dim(r)[1:2],type="mean",sPoly=as(q,"Spatial"))
  #mapMean <- mask(mapMean,q[q$NAME_1=="Québec",]) #spacePoly
  mapMean <- mask(mapMean,q) #spacePoly
  water<-predictors[["water"]]
  water[water<0.99]<-NA
  water<-resample(water,mapMean)
  mapMean <- mask(mapMean,water,inverse=TRUE) #spacePoly
  names(mapMean)<-paste(sp,names(m)[i])
  init[[ninit]]<-mapMean
  ninit<-ninit+1
  plot(mapMean, col = colo, axes = TRUE, box = FALSE, main = paste(names(m)[i],sp,sep=" - "))#,xlim=c(-1100,-900),ylim=c(200,300))
  #points(pressp,pch=16,col=alpha("forestgreen",0.99),cex=0.2)
  points(occsp,pch=16,col=alpha("black",0.95),cex=0.6)
  plot(st_geometry(q),add=TRUE,border=gray(0,0.15))
  ### contour
  ID<-inla.stack.index(attributes(m[[i]])$Stack,tag="est")$data
  #contour(mapMean,levels=quantile(m[[i]]$summary.fitted.values[["mean"]][ID],0.05),add=TRUE,col="red")
}
par(mfrow=c(1,1))

}


#####################################
### visualize spatial fields

#xlim<-range(st_coordinates(occ)[,1])
#ylim<-range(st_coordinates(occ)[,2])
xlim<-range(domain[,1])
ylim<-range(domain[,2])

proj<-inla.mesh.projector(Mesh,xlim=xlim,ylim=ylim,dims=c(300,300))

mfield<-inla.mesh.project(projector=proj,field=mpp$summary.random[['i']][['mean']])
sdfield<-inla.mesh.project(projector=proj,field=mpp$summary.random[['i']][['sd']])

par(mfrow=c(1,2),mar=c(3,3,2,5))

image.plot(list(x=proj$x,y=proj$y,z=mfield),col=viridis(100),asp=1,main="Spatial field") 
axis(1)
axis(2)
plot(st_geometry(q),add=TRUE,border="white")
plot(swe,add=TRUE,border=gray(0,0.5))
plot(occs,pch=1,cex=3*m$summary.fitted.values[index.est,"mean"],col=gray(0.1,0.3),add=TRUE)
brks<-c(0.01,0.25,0.50,0.75,0.99)
legend("topleft",pch=1,pt.cex=3*brks,col=gray(0,0.3),legend=brks,bty="n",title="Probability of location\nbeing an actual fire",inset=c(0.02,0.05))

image.plot(list(x=proj$x,y=proj$y,z=sdfield),col=viridis(100),asp=1,main="sd of spatial field (logit scale)") 
axis(1)
axis(2)
plot(swe,add=TRUE,border=gray(0,0.5))
plot(occs,pch=1,cex=3*m$summary.fitted.values[index.est,"mean"],col=gray(0.1,0.3),add=TRUE)
proj<-

