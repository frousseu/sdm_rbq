
### This file gathers predictors and other necessities such as shapefile for region delineation

### The output of thiis file is 
# na.shp
# predictors.tif
# predictors_transformed.tif

library(scales)
library(data.table)
library(rmapshaper)
library(spex)
library(dplyr)
library(raster)
library(RandomFields)
library(sf)
library(mapSpecies)
library(terra)
library(fasterize)

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


####################################################################
### raster of predictors and extent of study region
ext<-ext(c(-3500,3500,-3000,4000)) # North America
rp<-rast(ext,resolution=c(2,2),crs=prj)
na<-st_crop(na,ext)
st_write(na,"/data/predictors_sdm/na.shp",append=FALSE)
bb<-st_bbox(st_transform(na,4326))
plot(st_geometry(na))
plot(ext,add=TRUE)


##################################################
### WorldClim
# The following is to download the data
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
k<-c("tmean","prec","tmax","trange")
names(wc)<-k
wc<-project(wc,rp)
wc[[1]]<-wc[[1]]/10
wc[[2]]<-wc[[2]]/1000
wc[[3]]<-wc[[3]]/10
wc[[4]]<-wc[[4]]/10
for(i in k){
  wc[[i]]<-focal(wc[[i]],w=matrix(1,21,21),fun=mean,na.rm=TRUE,na.policy="only")
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

water<-mask(lc$water,st_transform(na,st_crs(lc)))
water <- subst(water, NA, 0)
ocean<-mask(lc$water,st_transform(na,st_crs(lc)),inverse=TRUE)
ocean <- subst(ocean, NA, 0)
names(ocean)<-"ocean"

lc<-lc[[names(lc)!="water"]]
lc<-c(lc,water,ocean)

# terrain
tri<-rast("/data/predictors_sdm/earthenv/terrain/tri_1KMmd_GMTEDmd.tif")
tri<-crop(tri,ext(bb[c(1,3,2,4)]))
tri<-aggregate(tri,2)
names(tri)<-"truggedness"

# elevation
elev<-rast("/data/predictors_sdm/earthenv/terrain/elevation_1KMmd_GMTEDmd.tif")
elev<-crop(elev,ext(bb[c(1,3,2,4)]))
elev<-aggregate(elev,2)
names(elev)<-"elevation"

# stack ee
ee<-c(c(lc,tri)/100,elev)
ee<-project(ee,rp)

###################################################
### Latitude and longitude
lat<-c(rp,rp)
lat<-setValues(lat,xyFromCell(lat,1:ncell(lat))[,1:2])
names(lat)<-c("longitude","latitude")


###################################################
### Ecoregions
# https://www.epa.gov/eco-research/ecoregions-north-america
#predictors<-rast("/data/predictors_sdm/predictors.tif")
ecoregions<-st_read("/data/predictors_sdm/ecoregions/NA_CEC_Eco_Level3.shp")
#ecoregions<-st_buffer(ecoregions,0)
ecoregions<-ms_simplify(ecoregions,keep=0.02)
ecoregions<-st_transform(ecoregions,st_crs(na))
#ecoregions<-st_intersection(ecoregions,na)
ecoregions <- st_union(ecoregions[,c("NA_L1NAME","geometry")],by_feature=TRUE)
ecoregions<-aggregate(ecoregions[,"geometry"],list(ecoregions$NA_L1NAME), st_union)
#eco2<-st_intersection(ecoregions,na)
l<-split(ecoregions,ecoregions$Group.1)

#rp<-predictors[[1]]
rpp<-raster(rp)
res<-lapply(l,function(i){
  rast(fasterize(i,rpp,background=0))
})
res<-res[sapply(res,function(i){max(values(i))})>0]
eco<-rast(res)
names(eco)<-paste0("eco",gsub(" |-|,|/","",names(res)))


#r<-predictors[["mixed"]]
#plot(r)
#plot(st_geometry(ecoregions),add=TRUE)
#e<-extract(r,vect(ecoregions),fun=mean,na.rm=TRUE)
#ecoregions$tmean<-e[,"tmean"]
#plot(ecoregions["NA_L1NAME"],nbreaks=200)
#plot(r)
#ecoregions$ecoreg<-ecoregions$Group.1
#plot(ecoregions["ecoreg"],reset=FALSE)
#text(st_coordinates(st_centroid(ecoregions)),#label=ecoregions$ecoreg,cex=0.75)


#############################################
### Stack all
predictors<-c(wc,ee,lat,eco)
png("/data/sdm_rbq/graphics/predictors.png",width=20,height=15,res=500,units="in")
plot(predictors,maxnl=100)
dev.off()

writeRaster(predictors,"/data/predictors_sdm/predictors.tif",overwrite=TRUE)

############################################
### Transform predictors ###################
############################################
# Do it here once it's faster than doing it before running models
predictors<-rast("/data/predictors_sdm/predictors.tif")

### Data transformations
# predictors<-scale(predictors)
# the following is the equivalent (see terra help)
means <- global(predictors, "mean", na.rm=TRUE)
pre <- predictors - means[,1]
sds <- global(pre, "rms", na.rm=TRUE)
predictors <- pre / sds[,1]
# function to backtransform to orginal scale
backTransform<-function(x,var){
  m<-match(var,row.names(means))
  me<-means[m,"mean"]
  sd<-sds[m,"rms"]
  (x*sd)+me
}

### remove factors temporarily
necos<-grep("eco",names(predictors))
ecos<-predictors[[necos]]
predictors<-predictors[[-necos]]

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

### add interactions with longitude
int<-names(predictors)[grep("conifers|cultivated",names(predictors))]
int<-int[-grep("_",int)]
int<-lapply(int,function(i){
  ras<-predictors[[i]]*predictors[["longitude"]]
  names(ras)<-paste0(i,"_","longitude")
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
predictors<-c(predictors,polys)

### bring back factors
predictors<-c(predictors,ecos)

writeRaster(predictors,"/data/predictors_sdm/predictors_transformed.tif",overwrite=TRUE)



