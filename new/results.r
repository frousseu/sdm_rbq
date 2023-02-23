
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
#ebird<-project(ebird,sdms[[1]])
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


