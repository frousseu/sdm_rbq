
# find . -type f -name '*Calidris_himantopus*'

# this is the part of results (cor, I) that comes from the old visualizations.r file

library(tidyr)
library(foreach)
library(doParallel)
library(FRutils)
library(berryFunctions)
library(magick)
library(jsonlite)
library(exiftoolr)
library(colorspace)
library(rmapshaper)
library(smoothr)

#######################################################################
### Get max value #####################################################
#######################################################################

lsdms<-list.files("/data/sdm_rbq/rasters",pattern="_birds.tif",full=TRUE)
#sdms<-lapply(lsdms,rast)

plan(multisession,workers=20)
cores<-nbrOfWorkers() # get nbr of workers from the chosen plan
chunks <- split(lsdms, rep(1:cores, each=ceiling(length(lsdms)/cores))[1:length(lsdms)])

maxs<-future_lapply(chunks,function(i){
  m<-lapply(i,function(j){
    global(rast(j,lyrs=1),"max",na.rm=TRUE)
  })
  m<-do.call("rbind",m)
  m$species<-basename(i)
  m
})
plan(sequential)
maxs<-do.call("rbind",maxs)
maxs$species<-gsub("_"," ",gsub("_birds.tif","",maxs$species))
df$max<-maxs$max[match(df$species,maxs$species)]


#######################################################################
### Get ratio value ###################################################
#######################################################################

# Since we already have max intensity value and convex hull area, we only need the convex hull of the predicted range

sp<-"Aratinga nenday"
occs<-getobs(sp)
sdm<-rast(file.path("/data/sdm_rbq/rasters/",paste0(gsub(" ","_",sp),"_birds.tif")))[[1]]



r<-sdm
r<-r/global(r,"max",na.rm=TRUE)[1,1]
th<-0.05
r[r<th]<-NA
r[r>=th]<-1
v <- as.polygons(r)
v<-st_as_sf(v)
vh<-st_convex_hull(v)
oh<-st_convex_hull(st_union(occs))
as.numeric(st_area(vh)/st_area(oh))


lsdms<-list.files("/data/sdm_rbq/rasters",pattern="_birds.tif",full=TRUE)
#sdms<-lapply(lsdms,rast)

plan(multisession,workers=15)
#cores<-nbrOfWorkers() # get nbr of workers from the chosen plan
#chunks <- split(lsdms, rep(1:cores, each=ceiling(length(lsdms)/cores))[1:length(lsdms)])
hr<-future_lapply(1:nrow(df),function(i){
   rpath<-file.path("/data/sdm_rbq/rasters",paste0(gsub(" ","_",df$species[i]),"_birds.tif"))
   r<-rast(rpath)[[1]]
   if(is.na(df$max[i]) | is.na(df$hullarea[i])){
     NA
   }else{
     r<-r/df$max[i]
     th<-0.1
     r[r<th]<-NA
     r[r>=th]<-1
     v <- as.polygons(r)
     v<-st_as_sf(v)
     vh<-st_convex_hull(v)
     oh<-df$hullarea[i]
     as.numeric(st_area(vh)/oh)
   }
})
plan(sequential)
df$hullratio<-unlist(hr)


table(df$reach<0.85)
table(df$reach>0.85 & df$hullratio>4)
table(df$reach>0.85 & df$hullratio<4 & df$max>10000)


#######################################################################
### Overlap with eBird ################################################
#######################################################################

lsdms<-list.files("/data/sdm_rbq/rasters",pattern="_birds.tif",full=TRUE)
lsdms<-lsdms[rev(order(file.info(lsdms)$mtime))]

df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
#sp<-gsub(" ","_",df$species[which(df$reach>=0.00)])
sp<-gsub(" ","_",unique(df$species))
lsdms<-lsdms[gsub("_birds.tif","",basename(lsdms))%in%sp]

lsdms<-lsdms[1:min(length(lsdms),5000)]
sdms<-lapply(lsdms,rast,lyrs="mean")
names(sdms)<-gsub("_birds.tif","",basename(lsdms))
dims<-sapply(sdms,function(i){paste(dim(i)[1:2],collapse="_")})
keep<-names(rev(sort(table(dims))))[1]
sdms<-rast(sdms[names(sdms)[dims%in%keep]])

matches<-match(gsub("_"," ",names(sdms)),d$species)
ebirdnames<-d$ebird[matches]
remove<-names(sdms)[is.na(ebirdnames)]
sdms<-sdms[[!names(sdms)%in%remove]] # remove what is not in ebird abundance 
ebirdnames<-d$ebird[match(gsub("_"," ",names(sdms)),d$species)]
ebirdnames<-gsub(" ","_",ebirdnames)

ebirdpath<-list.files("/data/predictors_sdm/expert_maps/eBird/abundance",full=TRUE,pattern="_ebird2.tif")
ebirdpath<-ebirdpath[unlist(sapply(ebirdnames,function(i){grep(i,ebirdpath)}),use.names=FALSE)]

ebird<-rast(lapply(ebirdpath,rast))
names(ebird)<-gsub("_ebird2.tif","",basename(ebirdpath))
# remove ebird models with values of 0 in NA
glob<-global(ebird,"max",na.rm=TRUE)
ebird<-ebird[[!names(ebird)%in%rownames(glob[glob$max==0,,drop=FALSE])]]

sdms<-sdms[[names(ebird)]]
#ebird<-project(ebird,crs(na))
#ebird<-crop(ebird,na)
#plot(ebird)

cors<-sapply(seq_len(nlyr(sdms)),function(i){
  print(i)
  #vsdms<-values(sdms[[i]])[,1]
  #vebird<-values(ebird[[i]])[,1]
  #cor(vsdms,vebird,use="complete.obs",method="pearson")
  #plot(vsdms,vebird)
  cor(values(sdms[[i]])[,1],values(ebird[[i]])[,1],use="complete.obs",method="pearson")
  #nicheOverlap(raster(sdms[[i]]),raster(ebird[[i]]),stat="I")
})
cors

#n<-10
#r1<-rast(matrix(1:n^2,ncol=n))#sdms[[i]]
#r2<-rast(matrix(1:n^2,ncol=n)+runif(n^2,10,100))#ebird[[i]]
#r1[sample(1:ncell(r1),70)]<-NA
#r2[sample(1:ncell(r2),70)]<-NA
#cor(values(r1)[,1],values(r2)[,1],use="complete.obs",method="pearson")
#layerCor(c(r1,r2),"pearson",na.rm=TRUE)


### nicheOverlap
cl<-makeCluster(21)
registerDoParallel(cl)
no<-foreach(i=names(sdms),.packages=c("dismo","raster")) %dopar% {
  r1<-raster(file.path("/data/sdm_rbq/rasters",paste0(i,"_birds.tif")),layer=1)
  r2<-raster(file.path("/data/predictors_sdm/expert_maps/eBird/abundance",paste0(i,"_ebird2.tif")),layer=1)
  nicheOverlap(r1,r2,stat="I")
}
stopCluster(cl)
no<-unlist(no)


### adds latest correlation to df of results
res<-data.frame(
  species=gsub("_"," ",names(sdms)),
  pearson=cors,
  I=no
)
df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
m<-match(df$species,res$species)
df$pearson<-ifelse(!is.na(m),res$pearson[m],df$pearson)
df$I<-ifelse(!is.na(m),res$I[m],df$I)
write.csv(df,file="/data/sdm_rbq/graphics/mapSpeciesres.csv",row.names=FALSE,append=FALSE)

########################################################################
### Rejected or not ####################################################
########################################################################

library(rmapshaper)
library(smoothr)

sp<-"Spizella passerina"
occs<-getobs(sp)
sdm<-rast(file.path("/data/sdm_rbq/rasters/",paste0(gsub(" ","_",sp),"_birds.tif")))[[1]]

r<-sdm
r<-r/global(r,"max",na.rm=TRUE)[1,1]
th<-0.1
r[r<th]<-NA
r[r>=th]<-1
v <- as.polygons(r)
v<-st_as_sf(v)
#v<-ms_simplify(v,0.01)
v<-smooth(v,method="ksmooth",smoothness=10)
par(mar=c(0,0,0,0))
plot(st_geometry(na),col="grey90",border=NA)
plot(st_geometry(v),col=adjustcolor("red",0.3),border=NA,lwd=5,add=TRUE)
naplot(lwd=0.5)
plot(st_geometry(occs),pch=16,cex=0.75,col=adjustcolor("black",0.5),add=TRUE)
vh<-st_convex_hull(v)
oh<-st_convex_hull(st_union(occs))
plot(st_geometry(vh),border=adjustcolor("red",0.3),lwd=5,add=TRUE)
plot(st_geometry(oh),border=adjustcolor("black",0.3),lwd=5,,add=TRUE)
as.numeric(st_area(vh)/st_area(oh))


########################################################################
### check vif for specific model #######################################
########################################################################

library(car)
vs<-unique(gsub("[[:digit:]]","",vars_pool))
vs<-vs[!vs%in%c("sbias","fixed")]
#vs<-c(vs,"water")
#vs<-vs[!vs%in%c("water_esa")]
dat<-as.data.frame(dmeshPred)
dat$y<-rnorm(nrow(dat))
mo<-lm(formula(paste("y~",paste(vs,collapse="+"))),data=dat)
vif(mo)

varset<-grep("_esa",vs,value=TRUE)
hist(rowSums(do.call("cbind",lapply(varset,function(i){
  backScale(dmeshPred[,i],i)
}))))
rowSums(lapply(dmeshPred[,grep("_esa",colnames(dmeshPred))]))


