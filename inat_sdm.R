
library(scales)
library(data.table)
library(rmapshaper)
library(spex)
library(dplyr)
library(raster)
library(RandomFields)
library(sf)
library(mapSpecies)
library(concaveman)
library(fields)
library(foreach) 
library(doParallel)
library(FRutils) 
library(viridisLite)
library(terra)
library(mapsf)

# set_ebirdst_access_key("vevof13s755j")


options(vsc.dev.args = list(width = 1200, height = 800))

inla.mesh2sp <- function(mesh) {
  crs <- inla.CRS(inla.CRSargs(mesh$crs))
  isgeocentric <- identical(inla.as.list.CRS(crs)[["proj"]], "geocent")
  if (isgeocentric || (mesh$manifold == "S2")) {
    stop(paste0(
      "'sp' doesn't support storing polygons in geocentric coordinates.\n",
      "Convert to a map projection with inla.spTransform() before
calling inla.mesh2sp()."))
  }
  
  triangles <- SpatialPolygonsDataFrame(
    Sr = SpatialPolygons(lapply(
      1:nrow(mesh$graph$tv),
      function(x) {
        tv <- mesh$graph$tv[x, , drop = TRUE]
        Polygons(list(Polygon(mesh$loc[tv[c(1, 3, 2, 1)],
                                       1:2,
                                       drop = FALSE])),
                 ID = x)
      }
    ),
    proj4string = crs
    ),
    data = as.data.frame(mesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
    match.ID = FALSE
  )
  vertices <- SpatialPoints(mesh$loc[, 1:2, drop = FALSE], proj4string = crs)
  
  list(triangles = triangles, vertices = vertices)
}

unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


#d<-fread("C:/Users/rouf1703/Documents/UdeS/Programmation/mapSpecies/observations-96473_20200626QC.csv")
#spinat<-d[,.N,by=.(taxon_species_name)][order(-N),][N>100,]$taxon_species_name
#d<-d[d$taxon_species_name%in%spinat,]
#d[,.(n=.N),by=.(taxon_species_name)][order(-n)][1:100]

#load("/data/predictors_sdm/predictors.RData")

prj<-"+proj=lcc +lat_0=49 +lon_0=-100 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
#predictors<-stack("/data/predictors_sdm/predictors.tif")
na<-st_read("/data/predictors_sdm/na.shp")

### load predictors
predictors<-rast("/data/predictors_sdm/predictors.tif")

### Data transformations
# predictors<-scale(predictors)
# the following is the equivalent (see terra help)
means <- global(predictors, "mean", na.rm=TRUE)
pre <- predictors - means[,1]
sds <- global(pre, "rms", na.rm=TRUE)
predictors <- pre / sds[,1]

### add quadratic terms
npredictors<-names(predictors)
predictors<-c(predictors,predictors^2)
names(predictors)<-c(npredictors,paste0(npredictors,2))
predictors<-predictors[[order(names(predictors))]]

### add interactions with latitude
int<-names(predictors)[!names(predictors)%in%c("latitude","latitude2","prec","prec2","tmean","tmean2","sbias")]
int<-lapply(int,function(i){
  ras<-predictors[[i]]*predictors[["latitude"]]
  names(ras)<-paste0(i,"_","latitude")
  ras
})
predictors<-c(predictors,rast(int))

### add interactions with tmean
int<-names(predictors)[!names(predictors)%in%c("latitude","latitude2","longitude","longitude2","tmean","tmean2","sbias")]
int<-int[-grep("_",int)]
int<-lapply(int,function(i){
  ras<-predictors[[i]]*predictors[[c("tmean","tmean2")]]
  names(ras)<-paste0(i,"_",c("tmean","tmean2"))
  ras
})
predictors<-c(predictors,rast(int))

### add latlon interactions
int<-names(predictors)[names(predictors)%in%c("latitude","latitude2")]
int<-lapply(int,function(i){
  ras<-predictors[["longitude"]]*predictors[[i]]
  names(ras)<-paste0("longitude","_",i)
  ras
})
predictors<-c(predictors,rast(int))


### add polynomials for tmean (test)
polys<-predictors[["tmean"]]^(3:6)
names(polys)<-paste0("tmean",3:6)
predictors<-c(predictors,rast(polys))
#predictors<-extend(predictors,ext(vect(st_buffer(na,1000))))
#for(i in names(predictors)){
#  predictors[[i]]<-focal(predictors[[i]],w=matrix(1,21,21),fun=mean,na.rm=TRUE,na.policy="only")
#  print(i)
#}
png("/data/sdm_rbq/graphics/extended_predictors.png",width=12,height=10,res=200,units="in")
plot(predictors[[2]])
plot(st_geometry(na),add=TRUE)
plot(Mesh,add=TRUE,vertex.color="grey50",edge.color=gray(0,0.1))
dev.off()


keep<-c("kingdom","phylum","class","order","family","genus","species","countryCode","stateProvince","decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters","day","month","year","recordedBy")
#d<-fread("D:/iNatGBIF/0137377-200613084148143.csv",encoding="UTF-8",select=keep)
#d<-fread("/data/predictors_sdm/inat/0137377-200613084148143.csv",encoding="UTF-8",select=keep) 
d<-fread("/data/predictors_sdm/inat/0212584-210914110416597.csv",encoding="UTF-8",select=keep) 
d<-d[countryCode%in%c("US","CA"),]
#prov<-c("Québec","Ontario","New Brunswick","Newfoundland and Labrador","Vermont","Maine","New Hampshire","New York","Massachusetts","Nunavut","Nova Scotia")
#prov<-c("Vermont","New Hampshire","Maine","Ontario","New York","Massachusetts","Québec","Rhode Island","Connecticut","New Brunswick","Nova Scotia")
prov<-c("Ontario","Québec","New Brunswick")
#d<-d[d$stateProvince%in%prov,]
#na<-na[na$NAME_1%in%prov,]
#d<-d[d$kingdom%in%c("Plantae"),] 
#d<-d[d$phylum%in%c("Chordata"),]
#d<-d[d$phylum%in%c("Arthropoda"),]
#d<-d[d$phylum%in%c("Tracheophyta"),]
d<-d[class%in%c("Aves"),]
#d<-d[d$class%in%c("Insecta"),]
#d<-d[d$order%in%c("Odonata"),]
#d<-d[d$class%in%c("Reptilia","Amphibia"),]
### get expert mps
#em1<-st_read("/data/predictors_sdm/expert_maps/IUCN/REPTILES.shp")
#em2<-st_read("/data/predictors_sdm/expert_maps/IUCN/AMPHIBIANS.shp")
#em<-st_read("/data/predictors_sdm/expert_maps/IUCN/FW_ODONATA.shp")
#em<-rbind(em1,em2)
#em$species<-em$binomial
#em<-st_read("/data/predictors_sdm/expert_maps/Little_USDA/Trees.shp")
em<-st_read("/data/predictors_sdm/expert_maps/eBird/eBird.shp")
em$species<-em$scntfc_
#em<-a[a$species%in%species,]
em<-st_transform(em,st_crs(na))


### subset dates
d[,date:=as.character(as.Date(paste(year,month,day,sep="-"),format="%Y-%m-%d"))]
d<-d[substr(date,6,10) >= "06-01" & substr(date,6,10) <= "08-01",]


d<-d[!is.na(decimalLatitude),]
d<-d[which(coordinateUncertaintyInMeters<50000),]
s<-st_as_sf(d,coords=c("decimalLongitude","decimalLatitude"),crs=4326)
s<-st_transform(s,st_crs(na))
#o<-st_intersects(s,na)

g<-st_make_grid(na,cellsize=10)
o<-st_intersects(s,g)
s$cell<-unlist(lapply(o,"[",1))
s<-s[!is.na(s$cell),]
s<-s[!duplicated(as.data.frame(s[,c("recordedBy","species","cell")])[,1:3]),]

#plot(st_geometry(na),border="grey85")
#plot(st_geometry(s),pch=16,col=gray(0,0.1),cex=0.5,add=TRUE)
#plot(log(r),axes=TRUE)
#plot(log(r),xlim=c(-400,800),ylim=c(-250,250))
#plot(st_geometry(na),border="grey65",add=TRUE)

e<-rasterize(vect(s[,1]),predictors[[1]],field=1,fun="length",background=0)
#e<-rasterize(s,predictors[[1]],field=1,fun="count",background=0)[[1]]
#e<-rasterize(s[s$genus%in%c("Carex"),1],predictors[[1]],field=1,fun="count",background=0)

sspecies<-s[,"species"]
cells<-raster::extract(e,vect(sspecies),cells=TRUE)
sspecies[,c("cell")]<-cells[,1]
sspecies<-sspecies[!duplicated(as.data.table(sspecies[,c("species","cell")])[,geometry:=NULL]),]
especies<-rasterize(vect(sspecies),predictors[[1]],field=1,fun="length",background=0)

#e2<-rasterize(occ,predictors[[1]],field=1,fun="count",background=0)
#raw<-rasterize(occ,rp[["sbias"]],field=1,fun="count",background=0)/rp[["sbias"]]

#save.image("/data/predictors_sdm/inat_sdm.RData")
#load("/data/predictors_sdm/inat_sdm.RData")

#################################################
#################################################

region<-st_buffer(concaveman(st_cast(na,"MULTIPOINT"),concavity=2),100)
#region<-st_buffer(concaveman(st_cast(s,"MULTIPOINT"),concavity=2),100)
region<-as(region,"Spatial")
region<-spTransform(region,crs(e)) ### added after for crs problem not sure if it works
#domain<-st_coordinates(st_cast(region,"MULTIPOINT"))[,1:2]
domain<-st_coordinates(st_cast(st_as_sf(region),"MULTIPOINT"))[,1:2]
domain<-st_coordinates(st_sample(st_as_sf(region),5000))
domain <- inla.nonconvex.hull(domain,convex = -0.015,resolution=75)

pedge<-0.01
edge<-min(c(diff(bbox(region)[0.2,])*pedge,diff(bbox(region)[2,])*pedge))
edge

png(file.path(getwd(),"graphics","mesh.png"),width=10,height=8,units="in",res=200)

#mesh<-inla.mesh.2d(loc.domain=NULL,max.edge=c(edge,3*edge),offset=c(edge,3*edge),cutoff=edge,boundary=domain,crs=CRS(proj4string(xs)))
Mesh <- inla.mesh.2d(loc.domain = NULL, #coordinates(occsp),
                     max.edge = c(edge,edge*3),
                     min.angle = 21,
                     cutoff = edge,
                     offset = c(edge,edge*3),
                     boundary=domain,
                     crs = crs(region))
plot(predictors[["tmean"]])
plot(Mesh,asp=1,add=TRUE)   
plot(region,add=TRUE)                  
dev.off()



#bpriors<-list(prec=list(default=1/(10)^2,barren=1/(10)^2,harsh=1/(10)^2,latitude=1/(10)^2,latitude2=1/(10)^2,Intercept=1/(100)^2,sbias=1/(100)^2),mean=list(default=0,Intercept=0,barren=0,harsh=0,latitude=0,latitude2=0,sbias=0))

bpriors<-list(prec=list(default=1/(1)^2,Intercept=1/(100)^2,sbias=1/(100)^2),mean=list(default=0,Intercept=0,sbias=0))

dummy<-setValues(predictors[[1:2]],runif(2*ncell(predictors[[1]])))
names(dummy)<-c("dummy1","dummy2")
predictors<-c(predictors,dummy)
rm(dummy)

#vars<-c("tmean","tmean2","latitude","latitude2","prec","builtup","cultivated","conifers","forested")
#vars<-c("latitude","latitude2","longitude","longitude2")
#vars<-c("forested","forested2","builtup","latitude","latitude2","tmean","tmean2")
#vars<-c("builtup","builtup2","cultivated","forested","forested2","cultivated2","prec","prec2","tmean","tmean2","herbaceous","shrubs")
#vars<-c("builtup","cultivated","forested","prec","tmean","tmean2","herbaceous","shrubs","broadleafs","mixed","conifers","water")
#vars<-c("tmean","tmean2","broadleafs","broadleafs_tmean","broadleafs_tmean2")
vars<-c("tmean","tmean2","tmean3","tmean4","tmean5")
#vars<-c("tmean","tmean2")
#vars<-c("1")
#vars<-c("tmax","tmax2","trange","builtup","forested")
#vars<-c("tmax","tmax2","forested","forested2","conifers","conifers2","builtup","builtup2","prec","prec2","truggedness","truggedness2")
#vars<-c("tmax","tmax2","forested","forested2","latitude","latitude2","barren","barren2","latitude:longitude","latitude2:longitude")
#vars<-c("tmax","tmax2","forested","forested2","barren")
#vars<-c("tmax","tmax2","forested","forested2","latitude","latitude2","barren")
#vars<-c("tmax","tmax2","forested","forested2","latitude","conifers","conifers2","builtup","harsh","truggedness")
#vars<-c("tmax","tmax2","conifers","conifers2","harsh","builtup")
#vars<-c("tmax","tmax2","forested","forested2","latitude")
#vars<-c("tmax","tmax2","latitude","latitude2","barren")
#vars<-c("tmax","tmax2")
#vars<-c("tmax")
#vars<-c("dummy1","dummy2")
#vars<-c("conifers","conifers2")
#vars<-c("latitude","latitude2") # problem with latitude????
#vars<-c("barren")
#vars<-names(predictors)
f<-formula(paste("y~",paste(vars,collapse="+")))



##############################################
##############################################

#species<-c("Bonasa umbellus","Archilochus colubris","Strix varia","Certhia americana","Melospiza lincolnii")
#species<-names(rev(sort(table(s$species)))[60:69])
species<-c("Setophaga fusca")
#species<-c("Acer saccharum","Acer saccharinum","Acer rubrum","Fraxinus nigra","Fraxinus americana","Fraxinus pennsylvanica","Populus balsamifera","Populus deltoides","Populus tremuloides","Populus grandidentata","Pinus banksiana","Pinus strobus","Picea mariana","Abies balsamea","Larix laricina")
#species<-c("Pinus banksiana","Picea mariana","Abies balsamea","Larix laricina","Tilia americana")
#species<-c("Picea mariana")
#species<-c("Agelaius phoeniceus","Cardinalis cardinalis","Turdus migratorius","Zenaida macroura")
#species<-c("Alces alces","Chelydra serpentina","Castor canadensis","Lepus americanus","Neovison vison")
#species<-c("Tamiasciurus hudsonicus")
#species<-names(rev(sort(table(s$species[s$phylum=="Tracheophyta"]))))[1:11]
#species<-names(rev(sort(table(s$species[s$class%in%c("Reptilia","Amphibia") & s$stateProvince=="Québec"]))))[1:20]
#species<-names(rev(sort(table(s$species[s$class%in%c("Reptilia") & s$stateProvince=="Québec"]))))[1:10]
#species<-names(rev(sort(table(s$species[s$phylum%in%c("Tracheophyta") & s$stateProvince=="Québec"]))))[1:10]
#species<-c("Arigomphus submedianus")
#species<-c("Picea mariana")
#species<-names(rev(sort(table(s$species[s$genus%in%c("Setophaga","Catharus") & s$stateProvince=="Québec"]))))[1:20]
#species<-names(rev(sort(table(s$species[s$class%in%c("Actinopterygii") & s$stateProvince=="Québec"]))))[1:20]
#species<-names(rev(sort(table(s$species[s$stateProvince=="Québec"]))))[1:21][-1]
#species<-names(rev(sort(table(s$species[s$class%in%c("Mammalia") & s$stateProvince=="Québec"]))))[1:20]
#species<-names(rev(sort(table(s$species))))[1:11]
#f<-function(i){as.integer(any(i==sp))}

init<-vector(mode="list",length=length(species)*2)
ninit<-1

loccs<-lapply(species,function(x){
  occ<-s[s$species==x,]  
  occsp<-as(occ,"Spatial")
})
names(loccs)<-species

predictorswrap<-wrap(predictors)
ewrap<-wrap(e)
especieswrap<-wrap(especies)

#rm(predictors,e,especies,dummy,dummy1,dummy2)

cl<-makeCluster(min(5,length(species)))
registerDoParallel(cl)
foreach(j=1:length(loccs),
        .packages=c("INLA","mapSpecies","raster","terra","sf","sp","exactextractr","scales","concaveman","viridis","mapsf"),
        .verbose=TRUE
        #.export=c("loccs","bpriors","f","na","predictorsw","em")
) %do% {
  
  #sink("/data/sdm_rbq/graphics/output.txt", append=TRUE)
  
  predictors<-rast(predictorswrap)
  e<-rast(ewrap)
  especies<-rast(especieswrap)
  
  #occs<-st_as_sf(loccs[[j]])
  #occs<-st_transform(occs,st_crs(r))
  
  sp<-names(loccs)[j] 
  occsp<-loccs[[j]]
  occs<-st_as_sf(occsp)
  
  eoccs<-rasterize(vect(occs),predictors[[1]],field=1,fun="length",background=0)
  eoccs<-setValues(eoccs,as.integer(values(eoccs)>0))
  vals<-values(e)*scales::rescale(values(eoccs/especies),to=c(1,max(values(especies))^1))
  vals<-ifelse(is.nan(vals) | is.infinite(vals),0,vals)
  esp<-setValues(e,vals)
  
  bdist<-500 #500
  buff<-st_buffer(concaveman(st_as_sf(st_cast(st_union(occs),"MULTIPOINT"),concavity=2)),bdist)
  add<-mask(esp,vect(buff),inverse=TRUE)
  #add[values(add)[,1]==0]<-1
  add[add==0]<-1
  #add<-mask(add,vect(na[na$NAME_0=="Canada" | na$NAME_1=="Alaska",]))
  add[is.na(add)]<-0
  esp<-esp+add
  
  #png("/data/sdm_rbq/graphics/esp.png",width=12,height=10,res=200,units="in")
  #plot(esp)
  #plot(log(aggregate(r[["sbias"]],2)),col=rev(viridis(200)))
  #plot(st_geometry(na),add=TRUE)
  #plot(Mesh,vertex.color="grey50",edge.color=gray(0,0.1),add=TRUE)
  #plot(mesh,col=gray(0,0.05),border="white",add=TRUE)
  #dev.off()
  
  names(esp)<-"sbias"
  rp<-c(predictors,esp)#/10000)
  
  r<-rp[[names(rp)[names(rp)%in%c(vars,"sbias")]]]
  #r<-aggregate(r,3)
  r<-extend(r,ext(vect(st_buffer(na,1000))))
  
  mesh<-st_buffer(st_union(st_as_sf(inla.mesh2sp(Mesh)$triangles)),100)
  test<-TRUE
  while(test){
    for(i in names(r)){
      r[[i]]<-focal(r[[i]],w=matrix(1,11,11),fun=mean,na.rm=TRUE,na.policy="only")
      #print(i)
    }
    vs<-extract(r[[1]],vect(mesh))
    nas<-is.na(vs[,1]) # used to be column = 2 why?
    test<-any(nas)
    print(sum(nas))
  }
  
  
  #png("/data/sdm_rbq/graphics/sbias.png",width=12,height=10,res=200,units="in")
  #plot(r[["sbias"]])
  #plot(log(aggregate(r[["sbias"]],2)),col=rev(viridis(200)))
  #plot(st_geometry(na),add=TRUE)
  #plot(Mesh,vertex.color="grey50",edge.color=gray(0,0.1),add=TRUE)
  #plot(mesh,col=gray(0,0.05),border="white",add=TRUE)
  #dev.off()
  
  
  #png("/data/sdm_rbq/graphics/explana.png",width=12,height=10,res=200,units="in")
  #plot(r[[1]])
  #plot(region,add=TRUE)
  #plot(Mesh,add=TRUE)
  #dev.off()
  
  #png("/data/sdm_rbq/graphics/use_predictors.png",width=12,height=10,res=200,units="in")
  #plot(r,maxnl=100)
  #plot(region,,add=TRUE)
  #plot(Mesh,add=TRUE)
  #dev.off()
  
  r<-stack(r)
  fixed<-r[["sbias"]]
  names(fixed)<-"fixed"
  r<-stack(r,fixed)
  explana<-explanaMesh(sPoly=region,meshSpace=Mesh,meshTime=NULL,X=r)
  
  weight <- ppWeight(sPoly = region, mesh = Mesh)
  f2<-update(f,reformulate(c(".","fixed")))
  #f<-y~tmean+mixed+conifers
  mpp <- ppSpace(f, sPoints = occsp,
                 explanaMesh = explana,
                 ppWeight = weight,
                 prior.range = c(50,0.1),
                 prior.sigma = c(1,0.1),
                 smooth = 2,
                 num.threads=2,
                 many=TRUE,
                 control.inla = list(int.strategy = "eb"),
                 fix = NULL,
                 sboffset = "sbias",
                 control.fixed = bpriors,
                 control.compute=list(config=TRUE),
                 verbose = TRUE
  ) 
  
  gb<-colorRampPalette(c("grey95","lightsteelblue1"))(1)
  
  colo<-list(
    #mean=colorRampPalette(c("grey90","steelblue","steelblue2","gold2","tomato2","red4"))(200),
    #mean=colorRampPalette(c("grey90","tomato","darkred","black"))(200)[1:200],
    mean=colorRampPalette(c(gb,"red4"))(200)[1:200],
    sd=rev(magma(200))[1:175]
  )
  
  png(paste0("/data/sdm_rbq/plots/","birds_",gsub(" ","_",sp),".png"),units="in",width=10,height=8,res=200)
  m<-list(mpp=mpp)
  par(mfrow=n2mfrow(length(m)),mar=c(0,0,0,4))
  for(i in 1:length(m)){
    type<-"mean"
    cols<-colo[[type]]
    mapMean<-mapSpace(m[[i]],dims=c(1,2)*dim(r)[1:2],sPoly=as(na,"Spatial"))#[[type]]
    #mapMean<-mapSpace(m[[i]],dims=c(1,2)*dim(r)[1:2],type=type,sPoly=as(na,"Spatial"))
    mapMean<-rast(mapMean)
    mapMean <- mask(mapMean,vect(na)) #spacePoly
    buff<-st_as_sf(st_buffer(st_convex_hull(st_union(occs)),dist=1000))
    mapMean <- mask(mapMean,vect(buff)) #spacePoly
    water<-predictors[["water"]]
    water[water<0.99]<-NA
    water<-resample(water,mapMean[[type]])
    mapMean <- mask(mapMean,water,inverse=TRUE) #spacePoly
    #names(mapMean)<-paste(sp,names(m)[i])
    init[[ninit]]<-mapMean
    ninit<-ninit+1
    plot(buff,border=FALSE)
    plot(mapMean[[type]],add=TRUE, col = cols, axes = TRUE, box = FALSE, main = paste(names(m)[i],sp,sep=" - "))#,xlim=c(-1100,-900),ylim=c(200,300))
    #plot(st_geometry(em[em$species==species[j],]),col=gray(0,0.2),border=NA,add=TRUE)
    #plot(hatchedLayer(st_geometry(em[em$species==species[j],]), "left2right",mode="sfc",density=5),col=gray(0,0.2),border="black",add=TRUE)
    points(occsp,pch=16,col=adjustcolor("black",0.25),cex=0.5)
    plot(st_geometry(na),add=TRUE,border=gray(0,0.15))
    ### contour
    ID<-inla.stack.index(attributes(m[[i]])$Stack,tag="est")$data
    if(type=="mean"){
      lim<-rasterToContour(raster(mapMean[[type]]),levels=max(values(mapMean[[type]]),na.rm=TRUE)*0.20)
      #plot(lim,add=TRUE,col=alpha("black",0.65),lwd=1)
    }
    crs(mapMean)<-crs(occsp)
    mtext(side=3,line=-5,adj=0.9,text=paste(sp,""),font=2,cex=1.4,col=gray(0,0.3))
    plot(attributes(weight)$dmesh,add=TRUE,border=gray(1,0.25),lwd=0.2)
    mf_scale()
    writeRaster(mapMean, filename=paste0("/data/sdm_rbq/rasters/",gsub(" ","_",sp),"_birds.tif"), overwrite=TRUE)
  }
  par(mfrow=c(1,1))
  dev.off()
  
  print(j) 
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
plot(st_geometry(na),add=TRUE,border="white")
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
  
  
  h<-raster("D:/NFI_MODIS250m_2011_kNN/NFI_MODIS250m_2011_kNN_Species_Frax_nig_v1.tif")
h<-aggregate(h,20)


raw<-rasterize(occ,rp[["sbias"]],field=1,fun="count",background=0)/rp[["sbias"]]
raw[raw<0.5]<-NA
raw<-rasterToPolygons(raw,na.rm=TRUE)
points(coordinates(raw)[,1],coordinates(raw)[,2],pch=1,cex=1,col="white")
plot(Mesh,add=TRUE,col=gray(0,0.9))

install.packages("/data/atlasBE", repos = NULL, type = "source")

install.packages("RPostgres")





inla.mesh2sp <- function(mesh) {
  crs <- inla.CRS(inla.CRSargs(mesh$crs))
  isgeocentric <- identical(inla.as.list.CRS(crs)[["proj"]], "geocent")
  if (isgeocentric || (mesh$manifold == "S2")) {
    stop(paste0(
      "'sp' doesn't support storing polygons in geocentric coordinates.\n",
      "Convert to a map projection with inla.spTransform() before
calling inla.mesh2sp()."))
  }
  
  triangles <- SpatialPolygonsDataFrame(
    Sr = SpatialPolygons(lapply(
      1:nrow(mesh$graph$tv),
      function(x) {
        tv <- mesh$graph$tv[x, , drop = TRUE]
        Polygons(list(Polygon(mesh$loc[tv[c(1, 3, 2, 1)],
                                       1:2,
                                       drop = FALSE])),
                 ID = x)
      }
    ),
    proj4string = crs
    ),
    data = as.data.frame(mesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
    match.ID = FALSE
  )
  vertices <- SpatialPoints(mesh$loc[, 1:2, drop = FALSE], proj4string = crs)
  
  list(triangles = triangles, vertices = vertices)
}



quebec<-readRDS("/data/predictors_sdm/worldclim/gadm36_CAN_1_sp.rds")
quebec<-quebec[quebec$NAME_1=="Québec",]
quebec<-st_as_sf(quebec)
quebec<-st_transform(quebec,prj)
quebec<-ms_simplify(quebec,keep=0.005)

mesh<-inla.mesh2sp(Mesh)
occs<-st_as_sf(loccs[[1]])
mesh<-st_as_sf(mesh$triangles)
r<-rp[["sbias"]]
occs<-st_transform(occs,st_crs(r))

e<-exact_extract(r,mesh)

mesh$effort<-sapply(e,function(i){
  sum(i[,1]*i[,2])	
})

plot(mesh["effort"])

o<-st_intersects(mesh,occs)

mesh$nobs<-lengths(o)
mesh$ratio<-mesh$nobs/mesh$effort
mesh$ratio<-ifelse(is.nan(mesh$ratio),NA,mesh$ratio)
mesh$cols<-colo.scale(mesh$ratio,cols=c("grey90","gold2","orange","red","red4"))
#mesh$cols<-colo.scale(sqrt(mesh$nobs),cols=c("grey90","gold2","orange","red","red4"))
mesh2<-st_crop(mesh,quebec)
#mesh2<-mesh

png(paste0("/data/sdm_rbq/plots/",gsub(" ","_",sp),"_effort.png"),units="in",width=10,height=10,res=300)
plot(st_geometry(mesh2),col=mesh2$cols,border=gray(0,0.2))
plot(st_geometry(quebec),add=TRUE,border="black",lwd=5)
plot(st_geometry(occs),add=TRUE,col=alpha("forestgreen",0.2),pch=16,cex=0.6)
dev.off()


#png(paste0("/data/sdm_rbq/plots/",gsub(" ","_",sp),"_effort.png"),units="in",width=10,height=10,res=300)
#plot(mesh2["effort"])
#plot(st_geometry(quebec),add=TRUE,border="black",lwd=5)
#plot(st_geometry(occs),add=TRUE,col=alpha("forestgreen",0.2),pch=16,cex=0.6)
#dev.off()

meshsf<-st_as_sf(mesh)
meshsf$tri<-1:nrow(meshsf)
w<-st_intersects(s,meshsf)
s$tri<-unlist(w)
trisp<-unique(as.data.table(s[,c("species","tri")])[,geometry:=NULL])
trisp<-trisp[,.(nsp=.N),by="tri"]
meshsf<-left_join(meshsf,trisp)
w<-st_intersects(meshsf,occs)
meshsf$pres<-as.integer(lengths(w)>0)
meshsf$prop<-meshsf$pres/meshsf$nsp
#meshsf2<-st_crop(meshsf,quebec)
meshsf2<-meshsf
meshsf2$cols<-colo.scale(meshsf2$prop,cols=c("grey90","gold2","orange","red","red4"),re=TRUE)

png(paste0("/data/sdm_rbq/plots/",gsub(" ","_",sp),"_prop.png"),units="in",width=10,height=10,res=300)
plot(st_geometry(meshsf2),col=meshsf2$cols)
plot(st_geometry(quebec),add=TRUE,border="black",lwd=5)
plot(st_geometry(occs),add=TRUE,col=alpha("forestgreen",0.2),pch=16,cex=0.6)
dev.off()



png(paste0("/data/sdm_rbq/plots/zzzz.png"),units="in",width=10,height=10,res=300,pointsize=10)
#plot(st_geometry(meshsf),lwd=0.1)
#plot(st_geometry(quebec),add=TRUE,border="black",lwd=1)
plot(aggregate(predictors[["harsh"]],13))
plot(st_geometry(meshsf),add=TRUE,lwd=0.1)
dev.off()




png(paste0("/data/sdm_rbq/plots/zzzz.png"),units="in",width=10,height=10,res=300)
#plot(e2)
plot(crop(e4,quebec))
plot(st_geometry(quebec),add=TRUE,border="black",lwd=1)
dev.off()


png(paste0("/data/sdm_rbq/plots/zzzz.png"),units="in",width=10,height=10,res=300)
h<-hist(values(esp),breaks=seq(0,16000,by=500),xlim=c(200,15000))
dev.off()





png(paste0("/data/sdm_rbq/plots/","zzzzzz.png"),units="in",width=24,height=12,res=200)
par(mfrow=c(2,4),mar=c(3,3,2,5))
xlim<-range(st_bbox(na)[c(1,3)])
ylim<-range(st_bbox(na)[c(2,4)])
proj<-inla.mesh.projector(Mesh,xlim=xlim,ylim=ylim,dims=c(300,300))
mfield<-inla.mesh.project(projector=proj,field=mpp$summary.random[['i']][['mean']])
sdfield<-inla.mesh.project(projector=proj,field=mpp$summary.random[['i']][['sd']])
image.plot(list(x=proj$x,y=proj$y,z=mfield),col=viridis(100),asp=1,main="Spatial field (log scale)") 
axis(1)
axis(2)
plot(st_geometry(na),add=TRUE,border=gray(0,0.5))
image.plot(list(x=proj$x,y=proj$y,z=sdfield),col=viridis(100),asp=1,main="sd of spatial field (log scale)") 
axis(1)
axis(2)
plot(st_geometry(na),add=TRUE,border=gray(0,0.5))
mean1<-rasterFromXYZ(cbind(as.matrix(expand.grid(x=proj$x,y=proj$y)),z=as.vector(mfield)),crs=crs(mapMean))
mean2<-resample(mean1,mapMean)
sd1<-rasterFromXYZ(cbind(as.matrix(expand.grid(x=proj$x,y=proj$y)),z=as.vector(sdfield)),crs=crs(mapMean))
sd2<-resample(sd1,mapMean)
plot(mapMean)
plot(exp(log(mapMean)-mean2))
plot(mask(mean2,na),col=colo.scale(seq(min(values(mean2),na.rm=TRUE),max(values(mean2),na.rm=TRUE),length.out=100),c("darkblue","blue","white","red","darkred"),center=TRUE))
plot(mask(sd2,na))
buff<-st_buffer(concaveman(st_as_sf(st_cast(st_union(occs),"MULTIPOINT"),concavity=2)),bdist)
plot(log(esp))
plot(st_geometry(na),add=TRUE)
plot(st_geometry(buff),add=TRUE)
plot(st_geometry(occs),cex=0.6,pch=1,col=alpha("black",0.5),add=TRUE)
add<-mask(esp,buff,inverse=TRUE)
add[values(add)==0]<-1
add<-mask(add,na[na$NAME_0=="Canada" | na$NAME_1=="Alaska",])
add[is.na(add)]<-0
plot(log(esp))
dev.off()


region<-st_buffer(concaveman(st_cast(na,"MULTIPOINT"),concavity=2),100)

dmesh<-st_as_sf(attributes(weight)$dmesh)
smean<-function(values,coverage){
  sum(values*coverage,na.rm=TRUE)/sum(coverage,na.rm=TRUE)
}
smeanbias<-function(values,coverage){
  sum(values*coverage,na.rm=TRUE)
}

r<-predictors[["builtup"]]
x<-exact_extract(r,dmesh,fun=smean,progress=TRUE)
dmesh$vals<-x

#r<-rp[["sbias"]]
r<-esp
x<-exact_extract(r,dmesh,fun=smeanbias,progress=TRUE)
dmesh$vals<-log(x)

plot(crop(r,vect(dmesh)))
plot(st_geometry(dmesh),add=TRUE,lwd=0.1)#,nbreaks=100,border=NA,key.pos=NULL)
plot(st_geometry(na),add=TRUE)
plot(dmesh,lwd=0.1,nbreaks=100,border=NA,key.pos=1,add=TRUE,pal=pal)

pal<-function(n){rev(magma(n,alpha=0.3))}
plot(st_buffer(dmesh,dist=-10)["vals"],lwd=0.1,nbreaks=100,border=NA,key.pos=NULL,pal=pal,add=TRUE)


reg<-na[na$NAME_1%in%c("Vermont","New Hampshire","Maine"),]
ps<-mask(mapMean,vect(reg))
ps<-crop(ps,vect(reg))
par(mar=c(0,0,0,0))
plot(ps)
plot(st_geometry(dmesh),add=TRUE,border="white",lwd=0.1)

plot(crop(predictors[["broadleafs"]],vect(reg)))
plot(crop(predictors[["broadleafs"]],vect(reg)))
plot(st_geometry(reg),add=TRUE)


x<-exact_extract(esp,dmesh,fun=smeanbias,progress=TRUE)
dmesh$esp<-log(x)
x<-exact_extract(esp+add,dmesh,fun=smeanbias,progress=TRUE)
dmesh$espadd<-log(x)
plot(dmesh[c("esp","espadd")])

plot(dmesh[c("esp")],key.pos=NULL)
plot(st_geometry(na),add=TRUE)


#####################################################
### add map pixelation

sd<-mapSpace(m[[i]],dims=c(1,1.5)*dim(r)[1:2],type="sd",sPoly=as(na,"Spatial"))

low<-mapSpace(m[[i]],dims=c(1,1.5)*dim(r)[1:2],type="0.025quant",sPoly=as(na,"Spatial"))

high<-mapSpace(m[[i]],dims=c(1,1.5)*dim(r)[1:2],type="0.975quant",sPoly=as(na,"Spatial"))

sdm<-c(rast(low),rast(high),rast(sd))
sdm<-log(sdm)
sdm<-mask(sdm,vect(buff))
sdm<-mask(sdm,vect(na))
crs(sdm)<-crs(predictors)
sdm<-aggregate(sdm,6)
#sdm<-disagg(sdm,2,method="bilinear")
v<-runif(ncell(sdm),values(sdm[[1]])[,1],values(sdm[[2]])[,1])
me<-values(mean(sdm[[1:2]]))
se<-values(diff(sdm[[1:2]])/4)
v<-rnorm(ncell(sdm),me,se)
table(is.na(v))
sdm2<-sdm[[1]]
sdm2<-setValues(sdm2,v)
#colos<-cols
colos<-colo.scale(200,c("white","lightgreen","green4","grey20"))
plot(exp(sdm2),col=colos)
plot(st_geometry(na),add=TRUE,border=gray(0,0.2))
plot(lim,add=TRUE,col=alpha("darkgreen",0.45),lwd=2)


# check CI with same color scale
par(mfrow=c(1,2),mar=c(0,0,0,0))
plot(sdm[[1]],range=range(values(sdm[[1:2]]),na.rm=TRUE))
plot(st_geometry(na),add=TRUE,border=gray(0,0.2))
plot(sdm[[2]],range=range(values(sdm[[1:2]]),na.rm=TRUE))
plot(st_geometry(na),add=TRUE,border=gray(0,0.2))



#################################################
### add range map predictor

#plot(st_geometry(occs),col="white")
#plot(predictors[["tmean"]])
#plot(st_geometry(na),add=TRUE)
#rm<-st_geometry(em[em$species==species[1],])
#plot(st_geometry(rm),col=gray(0,0.2),border=NA,add=TRUE)
#plot(st_geometry(occs),col=adjustcolor("firebrick",0.35),pch=16,add=TRUE)

#temp<-st_as_sf(crds(predictors[[1]],df=TRUE),coords=1:2)
#st_crs(temp)<-crs(predictors)
#dis<-st_distance(temp[1:100000,],rm)

test<-predictors[[1]]
test<-aggregate(test,5)
dis<-distance(test,vect(rm))
dis<-resample(dis,predictors[[1]])
dis<-setValues(dis,scales:::rescale(values(dis),0:1))
names(dis)<-"distrange"
predictors<-c(predictors,dis)

#predictors<-predictors[[-dim(predictors)[3]]


map_agreement<-function(x,y){
  if(!identical(st_crs(x),st_crs(y))){
    stop("x and y don't have the same crs")
  }
  #st_agr(x)<-"constant"
  #st_agr(y)<-"constant"
  int<-st_intersection(x,y)
  uni<-st_union(x,y)
  as.numeric(st_area(int)/st_area(uni))*ifelse(st_area(x)<=st_area(y),1,-1)
}


intensity2map<-function(x,pcutoff=0.20,pbandwidth=0.1){
  x[x<max(values(x),na.rm=TRUE)*pcutoff]<-NA
  x[!is.na(x)]<-1
  x<-as.polygons(x)
  #crs(map)<-prj
  x<-st_union(st_as_sf(x))
  #plot(st_geometry(x),border="grey90",add=TRUE)
  #x<-smooth(x,method="ksmooth",bandwidth=pbandwidth*mean(diff(st_bbox(x)[c(1,3,2,4)])[c(1,3)]))
  x
}

test<-intensity2map(mapMean[["mean"]],pcutoff=0.40,pbandwidth=0.01)

mesh<-inla.mesh2sp(Mesh)
#mesh<-attributes(weight)$mesh
dmesh<-attributes(weight)$dmesh


spem<-st_geometry(em[em$species==species[j],])
spem<-st_transform(spem,st_crs(test))
spem<-st_union(st_intersection(spem,na))

png("/data/sdm_rbq/graphics/mapOverlay.png",width=10,height=10,res=300,point=10,units="in")
plot(st_geometry(na))
plot(spem,border=NA,col=adjustcolor("blue",0.5),add=TRUE)
plot(test,border=NA,col=adjustcolor("red",0.5),add=TRUE)
legend("topright",legend=c("Expert map","Estimated map"),pch=15,cex=2,col=c(adjustcolor("blue",0.5),adjustcolor("red",0.5)),bty="n")
#plot(dmesh,add=TRUE)
dev.off()


map_agreement(st_buffer(test,0.0),spem)

png("/data/sdm_rbq/graphics/rangeHull.png",width=15,height=10,res=300,point=10,units="in")
rh<-rangeHull(occs,breaks=50)
dev.off()

png("/data/sdm_rbq/graphics/grid.png",width=10,height=10,res=300,point=10,units="in")
plot(st_geometry(na),col="grey90",border="white",lwd=2)
plot(dmesh,border=gray(0,0.1),add=TRUE)
plot(st_geometry(occs),pch=16,col="green4",add=TRUE)
dev.off()


smeanbias<-function(values,coverage){
  sum(values*coverage,na.rm=TRUE)
}

rb<-r[["sbias"]]
x<-exact_extract(rb,dmesh,fun=smeanbias,progress=TRUE)
dmesh$sbias<-x
plot(st_as_sf(dmesh))







class(mpp)<-"inla"
samps<-inla.posterior.sample(200,mpp)
samps<-lapply(samps,function(i){
  i$latent
})
samps<-do.call("cbind",samps)

unique(sapply(strsplit(row.names(samps),":"),"[",1))
table(sapply(strsplit(row.names(samps),":"),"[",1))

dims=c(1,2)*dim(r)[1:2]
type="mean"
sPoly=as(na,"Spatial")

mesh<-inla.mesh2sp(Mesh)$triangles

mapBasis <- inla.mesh.projector(attributes(mpp)$mesh,
                                dims = dims, xlim = c(xmin(sPoly), xmax(sPoly)),ylim = c(ymin(sPoly), ymax(sPoly)),crs = attributes(mpp)$mesh$crs)


ID <- inla.stack.index(attributes(mpp)$Stack, tag="pred")$data
mpp$summary.fitted.values[["mean"]][ID]
#samps[ID+4168,]

vals<-samps[grep("i:",row.names(samps)),]
vals<-samps[1:nrow(vals),]
#vals<-samps[(nrow(vals)+1):(nrow(vals)+1+nrow(vals)-1),]
vals<-rowMeans(vals)
#mapPred <- inla.mesh.project(mapBasis,modelSpace$summary.fitted.values[[type]][ID])
mapPred <- inla.mesh.project(mapBasis,vals)
mapRaster <- raster(t(mapPred[,ncol(mapPred):1]),xmn = min(mapBasis$x), xmx = max(mapBasis$x),ymn = min(mapBasis$y), ymx = max(mapBasis$y))

mapPred <- inla.mesh.project(mapBasis,cbind(vals,vals))
mat<-t(matrix(rev(mapPred[,1]),nrow=dims[1],ncol=dims[2],byrow=FALSE))
mapRaster <- raster(mat[,ncol(mat):1],xmn = min(mapBasis$x), xmx = max(mapBasis$x),ymn = min(mapBasis$y), ymx = max(mapBasis$y))
#mapRaster <- raster(t(as.vector(mapPred[,ncol(mapPred):1])),xmn = min(mapBasis$x), xmx = max(mapBasis$x),ymn = min(mapBasis$y), ymx = max(mapBasis$y))


plot(mapRaster,col=colo$mean)
plot(st_geometry(na),add=TRUE)
#plot(mesh,border=adjustcolor("white",0.5),lwd=0.2,add=TRUE)


preds<-mapSpace(m[[i]],dims=c(1,2)*dim(r)[1:2],sPoly=as(na,"Spatial"))
preds<-rast(preds)
preds<-mask(preds,vect(na))
dmesh<-st_as_sf(attributes(weight)$dmesh)

plot(exp(log(preds[["mean"]])-preds[["spacemean"]]),col=colo$mean)
plot(st_geometry(na),add=TRUE,lwd=0.5)
plot(st_geometry(dmesh),border=adjustcolor("black",0.05),add=TRUE)

plot(preds[["spacemean"]],col=colo$mean)


plot(mask(predictors[["truggedness"]],vect(na)))



dmesh<-st_as_sf(attributes(weight)$dmesh)
smean<-function(values,coverage){
  sum(values*coverage,na.rm=TRUE)/sum(coverage,na.rm=TRUE)
}

#r<-predictors[["builtup"]]
x<-exact_extract(predictors,dmesh,fun=smean,progress=TRUE)
dmesh$vals<-x

x<-exact_extract(predictors[[1:5]],dmesh)

### marginals

m<-mpp
class(m)<-"inla"
nsims<-2000
samples<-inla.posterior.sample(nsims,m)

params<-dimnames(m$model.matrix)[[2]]
nparams<-sapply(params,function(i){
  match(paste0(i,":1"),row.names(samples[[1]]$latent)) 
}) 
table(sapply(strsplit(row.names(samples[[1]]$latent),":"),"[",1))
nweights<-grep("i:",row.names(samples[[1]]$latent))

nd<-newdata2(as.data.frame(as.matrix(m$model.matrix)[,-1])[1:4120,],n=100,n2=3,fun=mean)
v<-names(nd)
nd<-unlist(nd,recursive=FALSE)

marginals<-lapply(seq_along(nd),function(k){
  p<-lapply(1:nsims,function(i){
    betas<-samples[[i]]$latent[nparams]
    names(betas)<-names(nparams)
    fixed<-cbind(Intercept=1,as.matrix(nd[[k]]))[,names(betas)] %*% betas # make sure betas and vars are in the same order
    # this if we want a spatial part
    #wk<-samples[[i]]$latent[nweights]
    #if(is.factor(xs@data[,v[k]])){ # factors never in model (et)
    #spatial<-as.matrix(inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,nlevels(size[,v[k]])),,drop=FALSE])) %*% wk
    #}else{
    #spatial<-as.matrix(AA) %*% wk # stack was Apn in fire
    #}
    #p<-fixed # ignores spatial part
    list(fixed,fixed)
  })
  p1<-do.call("cbind",lapply(p,"[[",1))
  p2<-do.call("cbind",lapply(p,"[[",2))
  p1<-t(apply(p1,1,function(i){c(quantile(i,0.025),mean(i),quantile(i,0.975))}))
  p2<-t(apply(p2,1,function(i){c(quantile(i,0.025),mean(i),quantile(i,0.975))}))
  p1<-exp(p1)
  p2<-exp(p2)
  p1<-setNames(as.data.frame(p1),c("low","mean","high"))
  add<-strsplit(names(nd)[k],"\\.")[[1]]
  data.frame(vars=add[1],block=add[2],p1)
})
names(marginals)<-names(nd)



preds<-cbind(do.call("rbind",marginals),do.call("rbind",nd))
preds<-split(preds,preds$vars)
preds<-lapply(preds,function(i){
  if(all(is.na(i$block))){
    list(i)
  }else{
    split(i,i$block)
  }
})


par(mfrow=n2mfrow(length(preds),asp=2),mar=c(4,4,1,1),oma=c(1,1,1,1))
lapply(names(preds),function(i){
  xvar<-strsplit(i,"_")[[1]][1]
  xvar2<-strsplit(i,"_")[[1]][2]
  xlim<-range(preds[[i]][[1]][,xvar])
  ywide<-unlist(preds)
  ylim<-range(as.numeric(ywide[grep("\\.mean",names(ywide))]))
  #ylim<-range(unlist(do.call("rbind",preds[[i]])[,c("low","mean","high"),drop=FALSE]))
  plot(0,0,xlab="",ylab="",xlim=xlim,ylim=ylim,type="n",mgp=c(2,0.5,0),tcl=-0.1)
  lapply(seq_along(preds[[i]]),function(jj){
    j<-preds[[i]][[jj]]
    x<-j[,xvar]
    lines(x,j[,"mean"],lwd=1.5,col=jj)
    polygon(c(x,rev(x),x[1]),c(j[,"low"],rev(j[,"high"]),j[,"low"][1]),col=alpha("black",0.1),border=NA)
    mtext(side=1,line=2,text=xvar,font=2,cex=1.25)
  })
  if(length(preds[[i]])>1){
    temp<-do.call("rbind",preds[[i]])
    var2<-strsplit(temp$vars[1],"_")[[1]][2]
    legend("top",title=var2,legend=round(as.numeric(unique(temp[,var2])),2),lwd=1,col=seq_along(preds[[i]]),bty="n",cex=1.25)
  }
  #p1<-marginals[[i]]
  #ylim<-c(0,max(unlist(marginals)))
  #ylim<-c(0,max(as.vector(p1[,2])))
  #plot(vals,p1[,2],type="l",ylim=ylim,xlab=i,ylab="Intensity",font=1,lty=1,mgp=c(1.5,0.45,0),tcl=-0.3)
  #lines(vals,p1[,2],lwd=1.5)
  #polygon(c(vals,rev(vals),vals[1]),c(p1[,1],rev(p1[,3]),p1[,1][1]),col=alpha("black",0.1),border=NA)
})

#preds1<-mapSpace(ml[[1]],dims=dim(r),sPoly=region)[["mean"]]



