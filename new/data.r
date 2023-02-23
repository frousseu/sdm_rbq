
### This file is responsibe for gathering the observations and the predictors 

source("/data/sdm_rbq/functions.r")

library(data.table)
library(terra)
library(sf)
library(jsonlite)
library(ebirdst)
library(tidyr)
library(berryFunctions)
library(FRutils)
library(rmapshaper)
#library(rnaturalearth)

options(width=150)
options(vsc.dev.args = list(width = 1000, height = 800))

#########################################################################
### Load predictors #####################################################

prj<-"+proj=lcc +lat_0=49 +lon_0=-100 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
#predictors<-stack("/data/predictors_sdm/predictors.tif")
na<-st_read("/data/predictors_sdm/na.shp")
nalakes<-st_read("/data/predictors_sdm/nalakes.gpkg")
nalakes<-st_transform(nalakes,st_crs(na))
areas<-as.numeric(st_area(nalakes))
ke<-paste(c("Manicouagan","-Jean"),collapse="|")
nalakes<-nalakes[unique(c(rev(order(areas))[1:50],grep(ke,nalakes$name_fr))),]
coast<-st_union(na)
na2<-ms_simplify(na,0.05)
nalakes2<-ms_simplify(nalakes,0.05)
coast2<-ms_simplify(coast,0.05)
#plot(st_geometry(na),col="grey90",border="white")
#plot(st_geometry(nalakes),col="lightblue",border=NA,add=TRUE)


### load predictors
#checkpoint("Predictor transformations")
op<-rast("/data/predictors_sdm/predictors.tif")

### Data transformations
# predictors<-scale(predictors)
# the following is the equivalent (see terra help)
#means <- global(op, "mean", na.rm=TRUE)
#pre <- op - means[,1]
#sds <- global(pre, "rms", na.rm=TRUE)
#predictors <- pre / sds[,1]
# function to backtransform to orginal scale
#backTransform<-function(x,var){
#  m<-match(var,row.names(means))
#  me<-means[m,"mean"]
#  sd<-sds[m,"rms"]
#  (x*sd)+me
#}

predictors<-rast("/data/predictors_sdm/predictors.tif")
#predictors<-rast("/data/predictors_sdm/predictors_transformed.tif")
#predictors<-aggregate(predictors,2)
#checkpoint("Predictor transformations done")

#######################################
### eBird info
ed<-as.data.table(ebirdst_runs)
ed<-ed[,.(common_name,scientific_name,species_code,resident,breeding_start,breeding_end)]
startend<-c("breeding_start","breeding_end")
ed[ ,(startend):=lapply(.SD,as.character),.SDcols=startend]
ed[,breeding_start:=ifelse(is.na(breeding_start),"2020-01-01",breeding_start)]
ed[,breeding_end:=ifelse(is.na(breeding_end),"2020-12-31",breeding_end)]
ed[,start:=substr(breeding_start,6,10)]
ed[,end:=substr(breeding_end,6,10)]
ed[,species:=scientific_name]
ed$species[ed$species=="Larus canus/brachyrhynchus"]<-"Larus brachyrhynchus"
#ed[,gbif:=ifelse(!is.na(match(species,unique(d$species))),species,NA)]

#k<-which(log(d$coordinateUn)>=10 & log(d$coordinateUn)<=10.5)
#hist(d$coordinateUn[k])
#plot(d$decimalLongitude[k],d$decimalLatitude[k])
#sample(d$occurrenceID[k],5)
#table(d$coordinatePrecision)

#################################################################
### Load data ###################################################

cols<-c("class","family","genus","species","infraspecificEpithet","countryCode","stateProvince","decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters","day","month","year","recordedBy","occurrenceID")
d<-fread("/data/predictors_sdm/inat/0425942-210914110416597.csv",encoding="UTF-8",nThread=10,select=cols) 
d<-d[class%in%c("Aves"),]
d<-d[countryCode%in%c("US","CA"),]
d<-d[!stateProvince%in%c("Hawaii"),]
d<-d[species!="",]


### Add family common names
fam<-unique(d$family)
fnames<-lapply(fam,function(i){
  print(i)
  x<-fromJSON(paste0("https://api.inaturalist.org/v1/search?q=",i,"&sources=taxa&per_page=30"))$results$record  
  x[which(x$matched_term==i),c("name","rank","matched_term","preferred_common_name")]
})
fnames<-do.call("rbind",fnames)
fnames<-setDT(cbind(family=fam,fnames))
fnames[,fname:=preferred_common_name]
fnames[,fname:=gsub(", and"," and",fname)]
d<-d[fnames[,c("family","fname")],on="family"]



### Manually replace certain names

# https://stackoverflow.com/questions/67908743/replace-column-values-in-table-with-values-from-lookup-based-on-matches-in-r-usi

#d$species[d$species=="Setophaga auduboni"]<-"Setophaga coronata"
#evening grosbeak
# larus canus
# barn owl
# double-crested cormorant
# american herring gull
# broad billed hummingbird

changes<-list(
  c("Setophaga auduboni","Setophaga coronata"),
  c("Phalacrocorax auritus","Nannopterum auritum"),
  c("Tyto furcata","Tyto alba"),
  c("Larus smithsonianus","Larus argentatus"),
  c("Larus brachyrhynchus","Larus canus/brachyrhynchus"),
  c("Grus canadensis","Antigone canadensis")
)
changes<-as.data.table(do.call("rbind",changes))
setnames(changes,c("old","new"))

d[changes, species := new, on = .(species = old)]


### Find iNaturalist names from obs IDs
d[,taxon:=trimws(paste(species,infraspecificEpithet))]
tax<-unique(d,by="taxon")[,.(taxon,occurrenceID)]
tax[,obsid:=basename(occurrenceID)]
ltax<-split(tax$obsid, ceiling(seq_len(nrow(tax))/30))


if(FALSE){ # not run if have been ran already
  taxanames<-lapply(seq_along(ltax),function(j){
    print(j)
    i<-ltax[[j]]
    x<-fromJSON(paste0("https://api.inaturalist.org/v1/observations/",paste(i,collapse=",")))
    taxaid<-x$results$taxon$id
    y<-fromJSON(paste0("https://api.inaturalist.org/v1/taxa?taxon_id=",paste(taxaid,collapse=","),"&per_page=500"))
    stopifnot(y$total_results<=500)
    res<-y$results[match(taxaid,y$results$id),c("name","preferred_common_name")]
    res<-res[match(i,x$results$id),]
    res
  })
  taxanames<-do.call("rbind",taxanames)
  names(taxanames)<-c("inat","inat_common")
  tax<-cbind(tax,taxanames)
  fwrite(tax,"tax.cvs",row.names=FALSE)
}

### check what happens to subspecies!!!!!! some are discarded
tax<-fread("tax.cvs")
#tax[,taxon:=gsub("Grus canadensis","Antigone canadensis",taxon)]
#fwrite(tax,"tax.cvs",row.names=FALSE)
d<-merge(d,tax[,.(taxon,inat,inat_common)])


#spnames<-d[,.(n=.N),by=.(species,taxon,inat,inat_common)]
#spnames[,matched:=inat%in%ed$species | species%in%ed$species | inat_common%in%ed$common_name | sub("^(\\S*\\s+\\S+).*", "\\1",inat)%in%ed$species]
#spnames<-spnames[order(matched,-n),]
#spnames[matched==FALSE,][1:100,]
#spnames[matched==TRUE,]

sprasters<-sapply(strsplit(list.files("/data/predictors_sdm/expert_maps/eBird/abundance"),"_"),function(i){paste(i[1:2],collapse=" ")})


### Get ebird names by matching with ebirdst dataset
spnames<-d[,.(n=.N),by=.(species,taxon,inat,inat_common)]
lnames<-lapply(unique(d$species),function(i){
  x<-spnames[species==i,.(species,taxon,inat,inat_common)]
  ma<-lapply(unlist(x),function(j){
    c(match(j,ed$species),match(j,ed$common_name))
  })
  ans<-unique(unlist(ma))
  ans<-ans[!is.na(ans)]
  if(length(ans)==0){
    NA
  }else{
    ed$species[ans]
  }
})
stopifnot(all(sapply(lnames,length)==1))
lnames<-data.table(species=unique(d$species),ebird=unlist(lnames))
d<-merge(d,lnames,by="species")


#############################
### keep screening

d<-d[!is.na(decimalLatitude),] # remove missing coordinates
d<-d[which(coordinateUncertaintyInMeters<50000),] # remove imprecise coordinates
d[,date:=as.character(as.Date(paste(year,month,day,sep="-"),format="%Y-%m-%d"))]
d[,md:=substr(date,6,10)]
d[,totobs:=.N,by=.(species)]
d<-d[totobs>10,]
d<-d[d$species!="",]
d[,.(n=.N),by=.(species,class)][n>200,][,.(n=.N),by=.(class)]
#prov<-c("Québec","Ontario","New Brunswick","Newfoundland and Labrador","Vermont","Maine","New Hampshire","New York","Massachusetts","Nunavut","Nova Scotia")
#prov<-c("Vermont","New Hampshire","Maine","Ontario","New York","Massachusetts","Québec","Rhode Island","Connecticut","New Brunswick","Nova Scotia")
#d<-d[d$stateProvince%in%prov,]

### get expert mps
#em<-st_read("/data/predictors_sdm/expert_maps/eBird/eBird.shp")
#em$species<-em$scntfc_
#em<-st_transform(em,st_crs(na))
#em1<-st_read("/data/predictors_sdm/expert_maps/IUCN/REPTILES.shp")
#em2<-st_read("/data/predictors_sdm/expert_maps/IUCN/AMPHIBIANS.shp")
#em<-st_read("/data/predictors_sdm/expert_maps/IUCN/FW_ODONATA.shp")
#em<-st_read("/data/predictors_sdm/expert_maps/Little_USDA/Trees.shp")


### subset dates
#d<-d[substr(date,6,10) >= "06-01" & substr(date,6,10) <= "07-25",]




### build spatial object
coords<-c("decimalLongitude","decimalLatitude")
s<-st_as_sf(d[,..coords],coords=coords,crs=4326)
s<-st_transform(s,st_crs(na))
d[,c("x","y"):=asplit(st_coordinates(s),2)]




### filter data
checkpoint("Make filter grid")
g<-rast(extent=ext(na)+5,resolution=5)
g<-setValues(g,1:ncell(g))
ex<-terra::extract(g,vect(s))
d$cell<-ex[,2]
d<-d[!is.na(cell),]
#d[,cell:=s$cell]
# slower older version
#system.time({
#  g<-st_make_grid(na,cellsize=5)
#  o<-st_intersects(s,g)
#  s$cell<-unlist(lapply(o,"[",1))
#})
checkpoint("Make filter grid done")

#s<-s[!duplicated(as.data.frame(s[,c("recordedBy","species","cell")])[,1:3]),]



#checkpoint("Rasterize for effort")
#e<-rasterize(vect(s[,1]),predictors[[1]],field=1,fun="length",background=0)
#sspecies<-s[,"species"]
#cells<-extract(e,vect(sspecies),cells=TRUE)
#sspecies[,c("cell")]<-cells[,"cell"]
#sspecies<-sspecies[!duplicated(as.data.table(sspecies[,c("species","cell")])[,geometry:=NULL]),]
#especies<-rasterize(vect(sspecies),predictors[[1]],field=1,fun="length",background=0)
#checkpoint("Rasterize for effort done")

