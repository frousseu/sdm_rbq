

library(sf)
library(ebirdst)
library(terra)
library(sf)
library(smoothr)
library(rnaturalearth)
library(data.table)
library(foreach)
library(doParallel)
library(rmapshaper)
library(gdalUtilities)

# set_ebirdst_access_key("1i058dhuuo4m",overwrite=TRUE)

cols<-c("class","family","genus","species","countryCode")
d<-fread("/data/predictors_sdm/inat/0212584-210914110416597.csv",encoding="UTF-8",select=cols) 
d<-d[class%in%c("Aves"),]
d<-d[countryCode%in%c("US","CA"),]
x<-d[,.(n=.N),by=.(species,family)][n>100,][rev(order(n)),]#[,.(n=.N),by=.(class)]


### North America
na<-st_read("/data/predictors_sdm/na.shp")
#proj <- "+proj=lcc +lat_0=49 +lon_0=-100 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
#na<-st_transform(na,proj)

ed<-as.data.table(ebirdst_runs)
ed<-ed[,.(common_name,scientific_name,species_code,resident,breeding_start,breeding_end)]
startend<-c("breeding_start","breeding_end")
ed[ ,(startend):=lapply(.SD,as.character),.SDcols=startend]

### missing species
writepath<-"/data/predictors_sdm/expert_maps/eBird/abundance"
already<-trimws(gsub("_ebird.tif|_"," ",list.files(writepath)))


#species_list<-x$species
species_list<-intersect(setdiff(d$ebird,already),ed$species)
species_list<-species_list[species_list!=""]
species_list<-species_list[sapply(gsub(" ","_",species_list),function(i){!any(grep(i,list.files(writepath)))})]

fileConn<-file("/data/sdm_rbq/graphics/output.txt")
cat("Hello",file="/data/sdm_rbq/graphics/output.txt",sep="\n")

species_list<-c("Calidris pusilla","Larus canus/brachyrhynchus")
#cl<-makeCluster(4)
#registerDoParallel(cl)
foreach(i=1:length(species_list),.packages=c("ebirdst","terra","sf")
) %do% {

  
  species<-species_list[i]
  spcode<-ed$species_code[ed$scientific_name==species]
  spcode

  sp_path <- ebirdst_download(species = spcode,tifs_only=TRUE,pattern="_abundance_seasonal_mean_hr_")
  abd <- load_raster(sp_path, product="abundance",period="seasonal",resolution="hr")
  w<-which(c("breeding","resident")%in%names(abd))
  if(length(w)>1){
    stop(paste(spcode,"without breeding or resident period"))  
  }
  abd<-rast(abd[[c("breeding","resident")[w]]])
  abd<-crop(abd,vect(st_transform(na,crs(abd))))
  abd<-mask(abd,vect(st_transform(na,crs(abd))))
  
  if(!dir.exists(writepath)){
    dir.create(writepath)
  }
  
  species<-gsub("canus/","",species)

  writeRaster(abd,filename=file.path(writepath,paste0(gsub(" ","_",species),"_ebird.tif")),overwrite=TRUE)

  folders<-list.files(tools::R_user_dir("ebirdst"),recursive=TRUE)
  g<-grep(spcode,folders)
  clean<-file.path(tools::R_user_dir("ebirdst"),folders[g])
  unlink(clean,recursive=TRUE)
  
  cat(paste(i,"of",length(species_list),"\n"),file="/data/sdm_rbq/graphics/output.txt", append=TRUE)

}
close(fileConn)
cat(paste(readLines("/data/sdm_rbq/graphics/output.txt"),"\n"))


list.files(tools::R_user_dir("ebirdst"),recursive=TRUE)
list.files(writepath)

if(FALSE){
  l<-list.files(writepath,full=TRUE)
  ebird<-rast(lapply(l,rast))
  names(ebird)<-gsub("_ebird.tif","",basename(l))

  #ebird<-project(ebird,crs(na))
  plot(ebird)
}


### Produce tifs with same res
# manually delete _ebird2.tif if any and recreate them
sdm<-rast("/data/sdm_rbq/rasters/Setophaga_citrina_birds.tif")
ebirds<-list.files(writepath,pattern="_ebird.tif",full=TRUE)
wkt<-st_crs(sdm)$wkt
te<-as.numeric(ext(sdm)[c(1,3,2,4)])
tr<-res(sdm)

if(any(grep("_ebird2.tif",list.files(writepath)))){
  stop("Remove \"_ebird2.tif\" first, otherwise gdalwarp does not work")
}
cl<-makeCluster(30)
registerDoParallel(cl)
foreach(i=seq_along(ebirds),.packages=c("gdalUtilities","sf","terra")) %dopar% {
  input <- ebirds[i]
  output <- gsub("_ebird.tif","_ebird2.tif",ebirds[i])
  gdalwarp(input,output,t_srs=wkt,te=te,tr=tr,r="bilinear")
}
stopCluster(cl)
