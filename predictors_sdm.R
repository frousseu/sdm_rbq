
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


colo<-colorRampPalette(c("grey90","steelblue4","steelblue2","gold","red1","red4"))(200)
prj<-"+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
#prj<-"+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +units=km +no_defs"

###################################################
### Québec 
q1<-raster:::getData('GADM', country='CAN', level=1)
q2<-raster:::getData('GADM', country='USA', level=1)
q1<-readRDS("/data/predictors_sdm/worldclim/gadm36_CAN_1_sp.rds")
q2<-readRDS("/data/predictors_sdm/worldclim/gadm36_USA_1_sp.rds")
q<-rbind(q1,q2)
q<-st_as_sf(q)
q<-st_transform(q,prj)
q<-ms_simplify(q,keep=0.005)


####################################################################
### raster of predictors and extent of study region
#ext<-extent(-2000,2000,-1000,3500)
#ext<-extent(as.vector(st_bbox(q)))
ext<-extent(-4500,2000,-2700,4600) # North America
rp<-raster(ext,resolution=c(5,5),crs=prj)
q<-st_crop(q,ext)
bb<-st_bbox(st_transform(q,4326))
plot(st_geometry(q))



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
#lon<-c(-105,-80,-55,-80,-105)
#lat<-c(45,45,45,70,70)
lon<-seq(bb[1],bb[3],length.out=5)
lat<-seq(bb[2],bb[4],length.out=3)
eg<-expand.grid(lon,lat)
lon<-eg[,1]
lat<-eg[,2]
l<-lapply(seq_along(lon),function(i){
  x<-raster::getData('worldclim',var='bio',res=0.5,lon=lon[i],lat=lat[i])
  x[[c(1,5,7,12)]]
})
#l<-list.files("predictors_sdm/wc0.5",pattern=".hdr|bio1|bio5|bio7|bio12")
#l<-list.files("predictors_sdm/wc0.5",pattern=".bil",full.names=TRUE)
#vs<-c("bio1","bio5","bio7","bio12")
#l<-lapply(l,raster)
WC<-do.call("merge",l)
wc<-WC
k<-c("tmean","tmax","trange","prec")
names(wc)<-k
wc<-projectRaster(wc,rp)
wc<-resample(wc,rp)
wc[[1]]<-wc[[1]]/10
wc[[2]]<-wc[[2]]/10
wc[[3]]<-wc[[3]]/10
wc[[4]]<-wc[[4]]/1000
for(i in k){
  wc[[i]]<-focal(wc[[i]],w=matrix(1,21,21),fun=mean,na.rm=TRUE,NAonly=TRUE)
}

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
lc<-stack(list.files(path,full.names=TRUE))
lc<-crop(lc,as.vector(bb)[c(1,3,2,4)])
lc<-aggregate(lc,2)
names(lc)<-names(habitats)[match(names(lc),paste0("Consensus_reduced_class_",unlist(habitats)))]
k<-names(habitats)
lc<-lc[[k]]
forested<-sum(lc[[c("conifers","broadleafs","mixed")]])
names(forested)<-"forested"
harsh<-sum(lc[[c("barren","ice")]])
names(harsh)<-"harsh"
lc<-stack(lc,forested,harsh)

# terrain
tri<-raster("C:/Users/rouf1703/Documents/UdeS/Programmation/mapSpecies/earthenv/terrain/tri_1KMmd_GMTEDmd.tif")
tri<-raster("/data/predictors_sdm/earthenv/terrain/tri_1KMmd_GMTEDmd.tif")
tri<-crop(tri,as.vector(bb)[c(1,3,2,4)])
tri<-aggregate(tri,2)
names(tri)<-"truggedness"

# stack ee
ee<-stack(lc,tri)/100
ee<-projectRaster(ee,rp)
ee<-resample(ee,rp)

# watermask
wm<-ee[["water"]]
wm[wm<0.999]<-NA
wm<-polygonize(wm)
watermask<-st_as_sf(st_union(wm))

###################################################
### Latitude

lat<-stack(rp,rp)
lat<-setValues(lat,coordinates(lat)[,1:2])
names(lat)<-c("longitude","latitude")

### Stack all
predictors<-stack(wc,ee,lat)
moytemp<-mean(values(wc$tmean),na.rm=TRUE)
sdtemp<-sd(values(wc$tmean),na.rm=TRUE)
predictors<-scale(predictors)
npredictors<-names(predictors)
predictors<-stack(predictors,predictors^2)
names(predictors)<-c(npredictors,paste0(npredictors,2))
predictors<-predictors[[order(names(predictors))]]
int<-names(predictors)[!names(predictors)%in%c("latitude","latitude2","prec","prec2","tmean","tmean2","sbias")]
int<-lapply(int,function(i){
  ras<-predictors[["latitude"]]*predictors[[i]]
  names(ras)<-paste0("latitude",":",i)
  ras
})
predictors<-stack(predictors,stack(int))
int<-names(predictors)[names(predictors)%in%c("latitude","latitude2")]
int<-lapply(int,function(i){
  ras<-predictors[["longitude"]]*predictors[[i]]
  names(ras)<-paste0("longitude",":",i)
  ras
})
predictors<-stack(predictors,stack(int))


##########################################################
### delete unnecessary objects
rm(list=ls()[!ls()%in%c("prj","colo","predictors","q","watermask")])
gc()

#save.image("/data/predictors_sdm/predictors.RData")

####################################################
####################################################
### viz

#png("/data/sdm_rbq/plots/zpredictors.png",units="in",width=10,height=8,res=75)
#plot(explana$X)#[[c("latitude")]])
#plot(mapMean)
#attributes(mpp)$X
#dev.off()
