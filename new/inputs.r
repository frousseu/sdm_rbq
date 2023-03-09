
library(mapSpecies)

### This where the list of species to submit is chosen along with data inputs to the model
### 
species<-c("Setophaga pinus","Setophaga americana","Setophaga petechia","Setophaga citrina","Setophaga ruticilla","Setophaga virens")
species<-c("Dolichonyx oryzivorus")
species<-c("Leucosticte australis","Toxostoma bendirei","Haematopus palliatus","Setophaga americana")#[4]
species<-c("Dryocopus pileatus","Setophaga coronata","Setophaga dominica")
tab<-table(d$species)
tab<-names(tab)[tab>30]
tab<-tab[d$ebird[match(tab,d$ebird)]%in%ed$scientific_name]
#tab<-tab[!d$ebird[match(tab,d$ebird)]%in%ed$scientific_name]
set.seed(sample(1:10000,9))
df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
#dontsp<-df$species[!is.na(df$I)]
dontsp<-df$species
tab<-tab[!tab%in%dontsp]
#species<-sample(tab,length(tab))
species<-species[sapply(species,function(i){print(i);nrow(getobs(i))>5})]
stopifnot(length(species)>0)
(lspecies<-species)
d[species%in%lspecies,][,.(n=.N),by=.(species)]

#loccs<-lapply(species,function(x){
#  occs<-s[s$species==x,]  
#  ### removes points with nn > 1000 km
#  nndist<-knn.dist(st_coordinates(occs),k=1)[,1]
#  occs<-occs[!nndist>=1200,]
#  occs
#})
#names(loccs)<-species

#plot(st_geometry(na))
#plot(st_geometry(loccs[[1]]),add=TRUE)

# https://becarioprecario.bitbucket.io/inla-gitbook/ch-priors.html
HN.prior = "expression:
  tau0 = 0.001;
  sigma = exp(-theta/2);
  log_dens = log(2) - 0.5 * log(2 * pi) + 0.5 * log(tau0);
  log_dens = log_dens - 0.5 * tau0 * sigma^2;
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);  
"
list(prec = list(prior = HN.prior))

bpriors<-list(prec=list(default=1/(0.5)^2,Intercept=1/(20)^2,sbias=1/(20)^2,fixed=1/(20)^2),mean=list(default=0,Intercept=0,sbias=0,fixed=0))

#colnames(dmeshPred)<-gsub(",","",colnames(dmeshPred))
#names(r)<-gsub(",","",names(r))
#names(predictors)<-gsub(",","",names(predictors))

###vars<-c("tmean","tmean2","broadleafs","builtup","cultivated","herbaceous","mixed","shrubs","conifers","harsh") #good

vars_pool<-c("tmean","tmean2","broadleafs","builtup","cultivated","herbaceous","mixed","shrubs","conifers","harsh") #good
#vars_pool<-c("tmean","tmean2","broadleafs","builtup","cultivated","herbaceous","forested","shrubs","conifers","harsh","water","logdistance")
vars_pool<-c("tmean","tmean2","broadleafs","builtup","cultivated","herbaceous","mixed","shrubs","conifers","harsh","water","logdistance")

vars_pool<-c("tmean","tmean2","deciduous_esa","builtup_esa","crop_esa","grass_esa","mixed_esa","shrubs_esa","conifers_esa","harsh_esa","logdistance","elevation","elevation2")

vars_pool<-c("tmean","tmean2","deciduous_esa","builtup_esa","crop_esa","grass_esa","mixed_esa","shrubs_esa","conifers_esa","logdistance","elevation","elevation2")

vars_pool<-c("tmean","tmean2","deciduous_esa","builtup_esa","crop_esa","grass_esa","mixed_esa","shrubs_esa","conifers_esa","logdistance","elevation")

vars_pool<-c("tmean","tmean2","deciduous_esa","builtup_esa","crop_esa","grass_esa","mixed_esa","shrubs_esa","conifers_esa","logdistance","elevation","elevation2")  # march one

vars_pool<-c("tmean","tmean2","deciduous_esa","builtup_esa","crop_esa","grass_esa","mixed_esa","shrubs_esa","conifers_esa","logdistance","elevation","elevation2","truggedness","water_esa")

#vars<-c("tmean","tmean2","broadleafs","builtup","cultivated","herbaceous","mixed","shrubs","conifers","harsh","conifers_longitude","longitude")
#vars<-c("forested","forested2")
#vars<-c("tmean","tmean2","broadleafs","builtup","cultivated","herbaceous","mixed","shrubs","conifers","harsh")
#vars<-c(grep("eco",colnames(dmeshPred),value=TRUE))
#vars<-colnames(dmeshPred)[-grep("eco",colnames(dmeshPred))]
#vars<-vars[1:42]
#vars<-vars[-grep("tmax|forested|barren|ice|flooded|mixed",vars)]
#vars<-c("tmean","tmean2","forested","forested2","elevation","elevation2","cultivated","conifers","longitude","conifers_longitude","cultivated_longitude")
#vars<-c("builtup")
#vars<-names(predictors)[-grep("3|4|5|6|tmax|tmin|trange|longitude|latitude|fixed|barren|flooded|ice|forested",names(predictors))]
#vars<-c("builtup")
#vars<-c("tmean","tmean2")

#os<-lapply(ls(),function(i){object.size(get(i))})
#os<-os[rev(order(unlist(os)))]
#lapply(os,format,units="Gb")

species

##########################################################
##########################################################
##########################################################
##########################################################

if(FALSE){

### From GBIF observations

# link to raw .tif found here
# https://io.biodiversite-quebec.ca/stac/collections/gbif_heatmaps/items/birds-heatmap

gbif<-rast("/data/predictors_sdm/gbif_heatmaps/gbif_birds_density_06-2022.tif")
gbif<-crop(gbif,vect(st_transform(dmesh,4326)))

target <- as.data.frame(rgbif::occ_data(scientificName = "Chordeiles minor",year="1900,2023",month="6,7",hasCoordinate = TRUE, occurrenceStatus="PRESENT",hasGeospatialIssue=FALSE, limit = 10000)$data)
#target<-target[target$month%in%6:7,]

target<-st_as_sf(target,coords=c("decimalLongitude","decimalLatitude"),crs=4326)
ov<-st_intersects(target,st_transform(dmesh,4326))
ov[sapply(ov,function(x) length(x)==0L)] <- NA
target$dmesh<-unlist(ov)
ex<-terra::extract(g,vect(st_transform(target,st_crs(na))))
target$cell<-ex[,2]
target<-target[!is.na(target$cell),]

plot(log(gbif))
plot(st_geometry(st_transform(na,4326)),add=TRUE)
plot(st_geometry(target),add=TRUE)


efftarget<-exact_extract(gbif,st_transform(dmesh,4326),fun="sum")
target<-st_transform(target,st_crs(na))
target$x<-st_coordinates(target)[,1] ###
target$y<-st_coordinates(target)[,2] ###
target<-as.data.table(target) ###
target<-unique(target,by=c("recordedBy","species","cell"))

}




