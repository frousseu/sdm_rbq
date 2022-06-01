
library(scales)
library(data.table)
library(rgbif)
library(rmapshaper)
library(spex)
library(dplyr)
library(raster)
library(RandomFields)
library(sf)
library(mapSpecies)
library(terra)

colo<-colorRampPalette(c("grey90","steelblue4","steelblue2","gold","red1","red4"))(200)
prj<-"+proj=lcc +lat_0=49 +lon_0=-100 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
#prj<-"+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +units=km +no_defs"

###################################################
### Québec 
can<-raster:::getData('GADM', country='CAN', level=1)
usa<-raster:::getData('GADM', country='USA', level=1)
can<-readRDS("/data/predictors_sdm/worldclim/gadm36_CAN_1_sp.rds")
usa<-readRDS("/data/predictors_sdm/worldclim/gadm36_USA_1_sp.rds")
na<-rbind(can,usa)
na<-st_as_sf(na)
na<-st_transform(na,prj)
na<-ms_simplify(na,keep=0.005)
st_write(na,"/data/predictors_sdm/na.shp")



####################################################################
### raster of predictors and extent of study region
#ext<-extent(-2000,2000,-1000,3500)
#ext<-extent(as.vector(st_bbox(na)))
ext<-ext(c(-3500,3500,-3000,4000)) # North America
#rp<-raster(ext,resolution=c(5,5),crs=prj)
rp<-rast(ext,resolution=c(2,2),crs=prj)
na<-st_crop(na,ext)
bb<-st_bbox(st_transform(na,4326))
plot(st_geometry(na))
plot(ext,add=TRUE)



####################################################
### Données OpenCanada
#h<-raster("C:/Users/rouf1703/Documents/UdeS/Programmation/mapSpecies/NFI_MODIS250m_2011_kNN_Species_Pice_Mar_v1.tif")
#h<-raster("C:/Users/rouf1703/Documents/UdeS/Programmation/mapSpecies/NFI_MODIS250m_2011_kNN_LandCover_VegNonTreed_v1.tif")
#h<-raster("C:/Users/rouf1703/Documents/UdeS/Programmation/mapSpecies/NFI_MODIS250m_2011_kNN_SpeciesGroups_Needleleaf_Spp_v1.tif")
#h<-raster("C:/Users/rouf1703/Documents/UdeS/Programmation/mapSpecies/NFI_MODIS250m_2011_kNN_SpeciesGroups_Broadleaf_Spp_v1.tif")
#ex<-extent(0.8e+06, 2.5e+06, 6.2e+06, 8.3e+06)
#h<-crop(h,ex)
#h<-aggregate(aggregate(aggregate(h,2),2),2)
#plot(h)


##################################################
### WorldClim
##lon<-c(-105,-80,-55,-80,-105)
##lat<-c(45,45,45,70,70)
#lon<-seq(bb[1],bb[3],length.out=5)
#lat<-seq(bb[2],bb[4],length.out=3)
#eg<-expand.grid(lon,lat)
#lon<-eg[,1]
#lat<-eg[,2]
#l<-lapply(seq_along(lon),function(i){
#  x<-raster::getData('worldclim',var='bio',res=0.5,lon=lon[i],lat=lat[i],path="../predictors_sdm/wc0.5")
#  x[[c(1,5,7,12)]]
#})

l<-list.files("../predictors_sdm/wc0.5",pattern=".bil",full.names=TRUE)
vs<-c("bio1_","bio5_","bio7_","bio12_")
l<-l[grep(paste(vs,collapse="|"),l)]
l<-split(l,sapply(strsplit(l,"_"),"[",2))
wc<-lapply(l,function(i){
  r<-lapply(i,function(j){
    rast(j)
  })
  do.call("merge",r)
})
wc<-rast(wc)
#l<-lapply(l,rast)
#WC<-do.call("merge",l[1:200])
#wc<-WC
k<-c("tmean","prec","tmax","trange")
names(wc)<-k
#wc<-projectRaster(wc,rp)
wc<-project(wc,rp)
#wc<-resample(wc,rast(rp))
wc[[1]]<-wc[[1]]/10
wc[[2]]<-wc[[2]]/1000
wc[[3]]<-wc[[3]]/10
wc[[4]]<-wc[[4]]/10
#for(i in k){
#  wc[[i]]<-focal(wc[[i]],w=matrix(1,21,21),fun=mean,na.rm=TRUE,NAonly=TRUE)
#}
for(i in k){
  wc[[i]]<-focal(wc[[i]],w=matrix(1,21,21),fun=mean,na.rm=TRUE,na.policy="only")
}
#png("/data/sdm_rbq/graphics/wc.png",width=12,height=10,res=200,units="in")
#plot(wc)
#dev.off()

###################################################
### Earth Env data
# https://data.earthenv.org/consensus_landcover/without_DISCover/Consensus_reduced_class_11.tif

path<-"C:/Users/rouf1703/Documents/UdeS/Programmation/mapSpecies/earthenv/landcover"
path<-"/data/predictors_sdm/earthenv/landcover"
list.files(path)

ext<-extent(-105,-50,40,84)

habitats<-list(
  conifers=1,
  broadleafs=3,
  mixed=4,
  shrubs=5,
  herbaceous=6,
  cultivated=7,
  flooded=8,
  builtup=9,
  ice=10,
  barren=11,
  water=12
)

# landcover
lc<-rast(list.files(path,full.names=TRUE))
lc<-crop(lc,ext(bb[c(1,3,2,4)]))
lc<-aggregate(lc,2)
names(lc)<-names(habitats)[match(names(lc),paste0("Consensus_reduced_class_",unlist(habitats)))]
k<-names(habitats)
lc<-lc[[k]]
forested<-sum(lc[[c("conifers","broadleafs","mixed")]])
names(forested)<-"forested"
harsh<-sum(lc[[c("barren","ice")]])
names(harsh)<-"harsh"
lc<-c(lc,forested,harsh)

# terrain
#tri<-rast("C:/Users/rouf1703/Documents/UdeS/Programmation/mapSpecies/earthenv/terrain/tri_1KMmd_GMTEDmd.tif")
tri<-rast("/data/predictors_sdm/earthenv/terrain/tri_1KMmd_GMTEDmd.tif")
tri<-crop(tri,ext(bb[c(1,3,2,4)]))
tri<-aggregate(tri,2)
names(tri)<-"truggedness"

# stack ee
ee<-c(lc,tri)/100
ee<-project(ee,rp)
#ee<-resample(ee,rp)

# watermask
wm<-ee[["water"]]
wm[wm<0.999]<-NA
wm<-polygonize(raster(wm))
watermask<-st_as_sf(st_union(wm))
st_write(watermask,"/data/predictors_sdm/watermask.shp")

###################################################
### Latitude

lat<-c(rp,rp)
#lat<-setValues(lat,1:ncell(lat))
lat<-setValues(lat,xyFromCell(lat,1:ncell(lat))[,1:2])
names(lat)<-c("longitude","latitude")

### Stack all
predictors<-c(wc,ee,lat)
png("/data/sdm_rbq/graphics/predictors.png",width=16,height=12,res=200,units="in")
plot(predictors,maxnl=100)
dev.off()
#moytemp<-mean(values(wc$tmean),na.rm=TRUE)
#sdtemp<-sd(values(wc$tmean),na.rm=TRUE)

### extend to cover future mesh
#k<-names(predictors)
#for(i in k){
#  predictors[[i]]<-focal(predictors[[i]],w=matrix(1,21,21),fun=mean,na.rm=TRUE,na.policy="only")
#}
#png("/data/sdm_rbq/graphics/wc.png",width=12,height=10,res=200,units="in")
#plot(wc)
#dev.off()

writeRaster(predictors,"/data/predictors_sdm/predictors.tif",overwrite=TRUE)



##########################################################
### delete unnecessary objects
#rm(list=ls()[!ls()%in%c("prj","colo","predictors","na","watermask")])
#gc()

#save.image("/data/predictors_sdm/predictors.RData")

####################################################
####################################################
### viz

#png("/data/sdm_rbq/plots/zpredictors.png",units="in",width=10,height=8,res=75)
#plot(explana$X)#[[c("latitude")]])
#plot(mapMean)
#attributes(mpp)$X
#dev.off()
