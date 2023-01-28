
library(INLA)
library(terra)
library(sf)
library(concaveman)
library(future)
library(future.apply)
library(FNN)
library(mapSpecies)
library(ebirdst)

### This file is used to build parameters and model which includes
### mesh, spde, priors, occ aggregation, model.matrix, stacks, etc.

### use more restricted zone
#bb<-st_transform(st_as_sf(st_as_sfc(st_bbox(c(xmin = -90, xmax = -57.6, ymax = 58, ymin = 42), crs = st_crs(4326)))),st_crs(predictors))
#bb<-st_intersection(bb,na)
#region<-st_buffer(concaveman(st_as_sf(st_sample(bb,5000)),concavity=2),50)
#plot(st_geometry(na))
#plot(st_geometry(region),border=NA,col=adjustcolor("black",0.1),add=TRUE)
#par(mar=c(0,0,0,0))
#plot(st_geometry(region))
#plot(crop(predictors[["elevation"]],ext(region)),add=TRUE)
#plot(st_geometry(na),add=TRUE)
#plot(st_geometry(occs),add=TRUE,cex=0.5)
#plot(st_geometry(region),add=TRUE)
#plot(Mesh,add=TRUE)


checkpoint("Region and domain")
region<-st_buffer(concaveman(st_cast(na,"MULTIPOINT"),concavity=2),50)
region<-st_transform(region,crs(predictors)) # was crs(e)### added after for crs problem not sure if it works
#domain<-st_coordinates(st_cast(region,"MULTIPOINT"))[,1:2]
#domain <- inla.nonconvex.hull(domain,convex = -0.015,resolution=75)
set.seed(1234)
domain<-st_coordinates(st_sample(region,5000))
domain <- inla.nonconvex.hull(domain,convex = -0.015,resolution=75)
checkpoint("Region and domain done")



pedge<-0.002
edge<-min(c(diff(st_bbox(region)[c(1,3)])*pedge,diff(st_bbox(region)[c(2,4)])*pedge))
edge


checkpoint("Mesh")
#mesh<-inla.mesh.2d(loc.domain=NULL,max.edge=c(edge,3*edge),offset=c(edge,3*edge),cutoff=edge,boundary=domain,crs=CRS(proj4string(xs)))
Mesh <- inla.mesh.2d(loc.domain = NULL, #coordinates(occsp),
                     max.edge = c(edge,edge*3),
                     min.angle = 21,
                     cutoff = edge,
                     offset = c(edge,edge*3),
                     boundary=domain,
                     crs = st_crs(na))
checkpoint("Mesh done")


png(file.path(getwd(),"graphics","mesh.png"),width=10,height=8,units="in",res=200)
plot(predictors[["tmean"]])
plot(Mesh,asp=1,add=TRUE)   
plot(st_geometry(region),add=TRUE)                  
dev.off()



checkpoint("Weights")
plan(multicore,workers=5)
weights <- ppWeight(sPoly = region, mesh = Mesh)
dmesh <- attributes(weights)$dmesh
dmesh$areas<-weights
dmesh$dmesh<-dmesh$id
checkpoint("Weights done")

plan(sequential)

#names(esp)<-"sbias"
#predictors<-predictors[[sapply(names(predictors),function(i){all(is.na(values(predictors[i])))})]]
r<-crop(predictors,vect(dmesh))#/10000)
r<-r[[!duplicated(names(r))]]
r<-r[[substr(names(r),1,3)!="eco"]]
#mesh<-st_convex_hull(st_buffer(st_union(st_as_sf(as.data.frame(Mesh$loc)[,1:2],coords=1:2)),50))
#r<-extend(r,ext(vect(st_buffer(mesh,50))))
#mesh<-st_convex_hull(st_buffer(st_union(st_as_sf(as.data.frame(Mesh$loc)[,1:2],coords=1:2)),50))
r<-extend(r,extend(ext(dmesh),res(r)))
crsr<-crs(r)

#########################################################
### Extend predictors to mesh ###########################
### Much faster nearest neighbour method for filling predictors
# https://gis.stackexchange.com/questions/433225/connecting-presence-point-to-nearest-raster-brick
# should modify code to only replace NA values and not all values
checkpoint("Nearest neighbour filling")
xy<-xyFromCell(r,1:ncell(r))
xy<-cbind(xy,values(r))
xynotna<-xy |> as.data.table() |> na.omit() |> as.matrix()
nn<-knnx.index(xynotna[,1:2],xy[,1:2],k=1)
r<-setValues(r,xynotna[nn,-(1:2)])
checkpoint("Done")
rm(xy)
rm(xynotna)
gc()
gc()


### Aggregate and focal method
# extends values with aggregated raster and then adds missing values to the same projection
#rr<-aggregate(r,50)
#test<-TRUE
#while(test){
#  rr<-focal(rr,w=matrix(1,3,3),fun=mean,na.rm=TRUE,na.policy="only")
#  vs<-terra::extract(rr[[1]],vect(mesh))
#  nas<-is.na(vs[,2]) # used to be column = 2 why?
#  test<-any(nas)
#  print(sum(nas))
#}
#rr<-resample(rr,r,method="bilinear")
#rr<-mask(rr,r,inverse=TRUE)
#r<-sum(rr,r,na.rm=TRUE)


#explana<-explanaMesh(sPoly=region,meshSpace=Mesh,meshTime=NULL,X=r)
#e1<-extract(r,Mesh$loc[,1:2])
#explana$Xmesh[,"tmean"]-e1[,"tmean"]




######################################################
### Extract predictors ###############################
### extract predictors to mesh
checkpoint("Extracting predictors for dual mesh cells")
plan(multicore,workers=10)
cores<-nbrOfWorkers() # get nbr of workers from the chosen plan
chunks <- split(1:nrow(dmesh), rep(1:cores, each=ceiling(nrow(dmesh)/cores))[1:nrow(dmesh)])
options(future.globals.maxSize = 1000 * 1024 ^ 2)
res<-future_lapply(chunks,function(chunksi){
t(exact_extract(r, 
                        dmesh[chunksi,], 
                        fun = function(values, coverage_fraction){
                        colSums(as.matrix(values) * coverage_fraction,
                                na.rm = TRUE) / sum(coverage_fraction)
                        },
                        force_df = FALSE,
                        progress = TRUE))
#res<-t(res)
})
plan(sequential)
dmeshPred <- do.call("rbind",res)
checkpoint("Done")


###############################################
### Add distance to coast

sea<-st_as_sfc(st_bbox(r))
sea<-st_difference(sea,st_union(na))
dis<-st_distance(sea,st_centroid(dmesh))
distance<-as.numeric(dis[1,])
logdistance<-log(distance+1,base=10)

#distance<-(distance-mean(distance))/sd(distance)
#adds<-matrix(mean(distance),ncol=1)
#colnames(adds)<-colnames(means)
#rownames(adds)<-"distance"
#means<-rbind(means,adds)
#adds<-matrix(sd(distance),ncol=1)
#colnames(adds)<-colnames(sds)
#rownames(adds)<-"distance"
#sds<-rbind(sds,adds)

#logdistance<-(logdistance-mean(logdistance))/sd(logdistance)
#adds<-matrix(mean(logdistance),ncol=1)
#colnames(adds)<-colnames(means)
#rownames(adds)<-"logdistance"
#means<-rbind(means,adds)
#adds<-matrix(sd(logdistance),ncol=1)
#colnames(adds)<-colnames(sds)
#rownames(adds)<-"logdistance"
#sds<-rbind(sds,adds)

dmeshPred<-cbind(dmeshPred,distance)
dmeshPred<-cbind(dmeshPred,logdistance)


#plot(dmesh["logdistance"],border=NA,nbreaks=200,reset=FALSE)
#plot(st_geometry(na),lwd=0.1,add=TRUE)
#plot(dmesh$distance,dmesh$logdistance)


######################################################
### Scale predictors  ################################
######################################################

#dmeshPred<-iris[,1:4]
means<-colMeans(dmeshPred)
sds<-apply(dmeshPred,2,sd)

Scale<-function(x){
  (x-mean(x))/sd(x)
}

backScale<-function(x,var){
  m<-match(var,names(means))
  me<-means[m]
  sd<-sds[m]
  (x*sd)+me
}

nps<-colnames(dmeshPred)
dmeshPred<-do.call("cbind",lapply(colnames(dmeshPred),function(i){
  Scale(dmeshPred[,i])
}))
colnames(dmeshPred)<-nps


######################################################
### Add predictor transformations ####################
######################################################

### add quadratic terms
trans<-dmeshPred^2
colnames(trans)<-paste0(colnames(trans),"2")
dmeshPred<-cbind(dmeshPred,trans)
dmeshPred<-dmeshPred[,order(colnames(dmeshPred))]


### add interactions with latitude
int<-colnames(dmeshPred)[!colnames(dmeshPred)%in%c("latitude","latitude2","prec","prec2","tmean","tmean2","sbias")]
l<-lapply(int,function(i){
  dmeshPred[,i]*dmeshPred[,"latitude"]
})
trans<-do.call("cbind",l)
colnames(trans)<-paste0(int,".","latitude")
dmeshPred<-cbind(dmeshPred,trans)

### add interactions with longitude
int<-colnames(dmeshPred)[grep("conifers|cultivated",colnames(dmeshPred))]
int<-int[-grep("\\.",int)]
l<-lapply(int,function(i){
  dmeshPred[,i]*dmeshPred[,"longitude"]
})
trans<-do.call("cbind",l)
colnames(trans)<-paste0(int,".","longitude")
dmeshPred<-cbind(dmeshPred,trans)

### add interactions with tmean
int<-colnames(dmeshPred)[!colnames(dmeshPred)%in%c("latitude","latitude2","longitude","longitude2","tmean","tmean2","sbias")]
int<-int[-grep("\\.",int)]
l<-lapply(int,function(i){
  res<-dmeshPred[,i]*dmeshPred[,c("tmean","tmean2")]
  colnames(res)<-paste0(i,".",c("tmean","tmean2"))
  res
})
trans<-do.call("cbind",l)
dmeshPred<-cbind(dmeshPred,trans)


### add latlon interactions
int<-colnames(dmeshPred)[colnames(dmeshPred)%in%c("latitude","latitude2")]
l<-lapply(int,function(i){
  res1<-dmeshPred[,"longitude"]*dmeshPred[,i]
  res2<-dmeshPred[,"longitude2"]*dmeshPred[,i]
  cbind(res1,res2)
})
trans<-do.call("cbind",l)
colnames(trans)<-c("longitude.latitude","longitude2.latitude","longitude.latitude2","longitude2.latitude2")
dmeshPred<-cbind(dmeshPred,trans)


### add polynomials for tmean (test)
l<-lapply(3:6,function(i){
  dmeshPred[,"tmean"]^i
})
trans<-do.call("cbind",l)
colnames(trans)<-paste0("tmean",3:6)
dmeshPred<-cbind(dmeshPred,trans)



######################################################
### Get dmesh ids
coords<-c("x","y")
s<-st_as_sf(d,coords=coords,crs=st_crs(na))
regionmesh<-region_mesh(Mesh) # modified region that uses the inner mesh boundary
o<-st_intersects(s,region) # keep only what's in the region (add region)
d<-d[lengths(o)==1L,] # keep only what's in the region (add region)
s<-st_as_sf(d,coords=coords,crs=st_crs(na)) # keep only what's in the region (add region)

o<-st_intersects(s,dmesh)
o[sapply(o,function(i){length(i)==0L})] <- NA
d[,dmesh:=unlist(o)]
rm(s);gc()

######################################################
### Recompute totobs
d[,totobs:=.N,by=.(species)]



#######################################################
### Effort info #######################################
#checkpoint("General info for effort")
#dmesh$nbobs<-lengths(st_intersects(dmesh,s))
#o<-unlist(st_intersects(s,dmesh),use.names=FALSE)
#dmeshdt<-as.data.table(dmesh)
#sdt<-as.data.table(s)
#sdt$id<-o
#sdt<-sdt[,.(nbspecies=length(unique(species))),by=.(id)]
#dmesh$nbsp<-merge(dmeshdt,sdt,all.x=TRUE)[is.na(nbspecies),nbspecies:=0][order(id),]$nbspecies
#checkpoint("Done")

#######################################################
### eBird periods #####################################

#https://stackoverflow.com/questions/41263946/remove-text-after-the-second-space


if(FALSE){

library(taxize)
library(rgbif)
library(ebirdst)

ebirdnames<-ed$species[is.na(ed$gbif)]
gbifnames<-unique(d$species)

l<-lapply(ebirdnames,function(i){
  print(i)
  x<-gnr_resolve(i)
  res<-x[x$data_source_title==c("GBIF Backbone Taxonomy"),]
  matched<-NULL
  if(nrow(res)>0){
    matched<-res$matched_name
  }
  x<-name_usage(name=i)$data
  if(nrow(x)>0){
    matched<-c(matched,x$species)
  }
  if(is.null(matched)){
    sp<-NA
  }else{
    matched<-sub("^(\\S*\\s+\\S+).*", "\\1",matched) #https://stackoverflow.com/questions/41263946/remove-text-after-the-second-space
    matched<-unique(matched)
    matched<-matched[!is.na(matched)]
    k<-matched%in%gbifnames
    if(any(k)){
      stopifnot(sum(k)==1L)
      sp<-matched[k]
    }else{
      x<-name_lookup(query=i)$data
      x<-x[x$class%in%"Aves",]
      matched<-x$species  
      matched<-sub("^(\\S*\\s+\\S+).*", "\\1",matched) #https://stackoverflow.com/questions/41263946/remove-text-after-the-second-space
      matched<-unique(matched)
      matched<-matched[!is.na(matched)]
      k<-matched%in%gbifnames
      if(any(k)){
        #stopifnot(sum(k)==1L)
        if(sum(k)>1){
          sp<-NA # paste(matched[k],collapse="/")
        }else{
          sp<-matched[k]
        }
      }else{
        sp<-NA
      }
    }
  }
  res<-c(species=i,species_match=sp)
  #print(res)
  res
})
x<-as.data.table(do.call("rbind",l))

k<-is.na(ed$gbif)
ed$gbif[k]<-x$species_match[match(ed$species[k],x$species)]


ebirdnames<-ed$species[is.na(ed$gbif)]

miss<-setdiff(gbifnames,ed$gbif)
d[species%in%miss,.(n=.N),by=species][order(-n),][1:50,]



sp<-c("Setophaga coronata","Setophaga auduboni")
sp<-c("Setophaga auduboni")
sp<-c("Setophaga coronata")
s<-st_as_sf(d[species%in%sp,],coords=c("x","y"),crs=st_crs(na))
plot(st_geometry(na))
plot(st_geometry(s),add=TRUE)




#writepath<-"/data/predictors_sdm/expert_maps/eBird/abundance"
#x<-list.files(writepath)


grep("/",ebirdnames)

##### find names
library(rgbif)
library(taxize)


#spnames<-unique(d$species[!d$species%in%ed$species])
l<-lapply(spnames,function(i){
  print(i)
  x <- gnr_resolve(i)
  res<-x[x$data_source_title==c("The eBird/Clements Checklist of Birds of the World"),]
  if(nrow(res)>0){
    matched<-sub("^(\\S*\\s+\\S+).*", "\\1",res$matched_name) #https://stackoverflow.com/questions/41263946/remove-text-after-the-second-space
    k<-matched%in%gbifnames
    if(any(k)){
      stopifnot(sum(k)==1L)
      sp<-matched[k]
    }else{
      key<-name_backbone(i)$usageKey
      x<-name_usage(key=key,data="synonyms")$data$canonicalName
      matched<-unique(name_usage(name=i)$data$species)
      matched<-matched[!is.na(matched)]
      stopifnot(length(matched)>=1L)
      matched<-sub("^(\\S*\\s+\\S+).*", "\\1",matched)
      k<-matched%in%gbifnames
      if(any(k)){
        stopifnot(sum(k)==1L)
        sp<-matched[k]
      }else{
        sp<-NA
      }
    }
    c(species=i,species_match=sp)
  }
})
x<-as.data.table(do.call("rbind",l))


as.data.table(name_lookup("Circus hudsonius"))[,..int]



spnames<-unique(d$species)

l<-lapply(ed$species[1:10],function(i){
  print(i)
  key<-name_backbone(i)$usageKey
  #print(key)
  x<-name_usage(key=key,data="synonyms")$data
  
#x[,c("scientificName","species","canonicalName","accepted","rank")]
#grep(" x ",as.data.table(name_lookup(query=s)$data)$species)

  if(nrow(x)==0){
    x<-as.data.table(name_lookup(query=i)$data)$species
    sp<-NA
  }else{
    sp<-x$canonicalName
    sp<-unique(sapply(strsplit(sp," "),function(j){paste(j[1:2],collapse=" ")}))
    k<-sp%in%spnames
    if(!any(k)){
      sp<-NA
    }else{
      stopifnot(sum(k)==1)
      sp<-sp[k]
    }
  }
  sp
})
 

spn<- c("scientificName","accepted","species","canonicalName")

i<-"Circus cyaneus"

x<-unique(as.data.table(name_lookup(query=i)$data)$canonicalName)

i<-"Circus hudsonius"
key<-name_backbone(i)$usageKey
x<-name_usage(key=key,data="synonyms")$data
x$species

i<-
x<-as.data.table(name_usage(name=i)$data)
x<-as.data.table(name_lookup(query=i)$data)

int<-intersect(spn,names(x))
x[,..int]
x$species

}


### Create result file
#df<-data.frame(matrix(ncol = 10, nrow = 0))
#colnames(df)<-c("species","date","n","reach","predictors","range","sd","pearson","spearman","I","D")
#write.table(df,file="/data/sdm_rbq/graphics/mapSpeciesres.csv",row.names=FALSE,header=TRUE,sep=",",append=FALSE)
#df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")


#############################################################
### Mapping predictors ######################################

# optional mapping of predictors with original predictors op

if(FALSE){
    checkpoint("Mapping predictors")
    mf<-n2mfrow(length(names(op)),asp=4)
    dmesh<-st_as_sf(attributes(weights)$dmesh)
    msmean<-function(values,coverage){
    colSums(values*coverage,na.rm=TRUE)/sum(coverage,na.rm=TRUE)
    }
    x<-t(exact_extract(op[[names(op)]],dmesh,fun=msmean,progress=TRUE,max_cells_in_memory=3e+12))
    dmesh[,names(op)]<-x
    png("/data/sdm_rbq/graphics/predictors_dmesh.png",width=mf[2]*4,height=mf[1]*4,res=500,units="in")
    par(mfrow=mf,bg="black")
    lapply(names(op),function(i){
    plot(dmesh[i],pal=magma,nbreaks=100,lwd=0.05,border=adjustcolor("white",0.25),mar=c(0,0,0,0),key.pos=NULL,reset=FALSE)
    plot(st_geometry(na),lwd=0.30,border=adjustcolor("white",0.85),mar=c(0,0,0,0),add=TRUE)
    mtext(side=3,line=-6,text=i,adj=0.90,col="white",font=2)
    })
    dev.off()
}


#options(vsc.dev.args = list(width = 1200, height = 600))