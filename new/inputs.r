
library(mapSpecies)

### This where the list of species to submit is chosen along with data inputs to the model
### 
species<-c("Artemisiospiza belli","Charadrius nivosus","Lagopus lagopus","Lanius borealis","Stercorarius parasiticus","Thalasseus sandvicensis","Tyrannus melancholicus")
tab<-table(d$species)
tab<-names(tab)[tab>30]
tab<-tab[d$ebird[match(tab,d$ebird)]%in%ed$scientific_name]
set.seed(sample(1:10000,1))
df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
tab<-tab[!tab%in%df$species]
#species<-sample(tab,9)
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


bpriors<-list(prec=list(default=1/(1)^2,Intercept=1/(20)^2,sbias=1/(20)^2,fixed=1/(20)^2),mean=list(default=0,Intercept=0,sbias=0,fixed=0))

#colnames(dmeshPred)<-gsub(",","",colnames(dmeshPred))
#names(r)<-gsub(",","",names(r))
#names(predictors)<-gsub(",","",names(predictors))

###vars<-c("tmean","tmean2","broadleafs","builtup","cultivated","herbaceous","mixed","shrubs","conifers","harsh") #good

vars_pool<-c("tmean","tmean2","broadleafs","builtup","cultivated","herbaceous","mixed","shrubs","conifers","harsh") #good
vars_pool<-c("tmean","tmean2","broadleafs","builtup","cultivated","herbaceous","mixed","shrubs","conifers","harsh","water","logdistance")
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

