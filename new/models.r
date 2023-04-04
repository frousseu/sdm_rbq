 

# nohup bash run.sh > nohup.out 2>&1 &

source("/data/sdm_rbq/functions.r")

if(FALSE){
  source("/data/sdm_rbq/functions.r")
  source("/data/sdm_rbq/data.r")
  source("/data/sdm_rbq/parameters.r")
  source("/data/sdm_rbq/inputs.r",echo=TRUE)
}

library(terra)
library(sf)
library(data.table)
library(concaveman)
library(future)
library(future.apply)
library(future.callr)
library(viridisLite)
library(magick)
library(FNN)
library(INLA)
library(mapSpecies)
library(progressr)

print(paste("Loading data",Sys.time(),collapse=" "))

load("/data/sdm_rbq/parameters.RData")
op<-rast(opwrap)
predictors<-rast(predictorswrap)
r<-rast(rwrap)
source("/data/sdm_rbq/functions.r",echo=TRUE)
source("/data/sdm_rbq/inputs.r",echo=TRUE)

checkpoint("Data loaded")
options(vsc.dev.args = list(width = 1000, height = 800))

### This runs models for each species and produces necessary model outputs
# The number of cores used is adjusted according to the number of species

inla.setOption(pardiso.license = "/home/rouf1703/pardiso.lic")

dims<-dim(r)[1:2]

#cl<-makeCluster(min(2,length(species)))
#registerDoParallel(cl)

#predictors<-rast(predictorswrap)
#e<-rast(ewrap)
#especies<-rast(especieswrap)

#foreach(j=1:length(loccs),
#  .packages=c("INLA","mapSpecies","raster","terra","sf","sp","exactextractr","scales","concaveman","viridis","mapsf"),
#  .verbose=TRUE
  #.export=c("loccs","bpriors","f","na","predictorsw","em")
#) %dopar% {

#plan(sequential)
plan(multisession,workers=min(c(3,length(species))))
cores<-nbrOfWorkers() # get nbr of workers from the chosen plan
#plan(sequential)
options(future.globals.maxSize = 5000 * 1024 ^ 2)
options(progressr.enable=TRUE)

print(paste(species,collapse=" "))

ml<-future_lapply(seq_along(species),function(j){

#j<-1
#prog<-progressr::progressor(along = species) # progressr messages not working
#prog(message = sprintf("Starting %s", species[j]))

sp<-species[j] 
print(paste(Sys.time(),"Running species",sp,collapse=" "))
be<-unlist(ed[match(sp,ed$scientific_name),c("start","end")])

### all obs for period
if(be[1]<be[2]){
  obs<-d[md>=be[1] & md<=be[2],]
}else{
  # Tyrannus savanna breeds from nov to jan
  dr<-substr(seq.Date(as.Date(paste0("2007-",be[1])),as.Date(paste0("2008-",be[2])),by=1),6,10)
  obs<-d[md%in%dr,]
}

#obs<-d

### remove duplicate obs
obs<-unique(obs,by=c("recordedBy","species","cell"))
### remove places
#sb<-st_as_sf(obs,coords=c("x","y"),crs=st_crs(na)) 
#osb<-st_intersects(sb,na[na$NAME_1%in%c("Ontario","Québec"),])
#osb<-st_intersects(sb,na[na$NAME_0%in%c("Canada"),])
#obs<-obs[as.logical(lengths(osb)),]
#writeRaster(aggregate(sdm[[c("mean","linksd")]],2), filename=paste0("/data/sdm_rbq/Chordeiles_minor.tif"), overwrite=TRUE)
### keep certain n of obs
#w<-which(obs$species==sp)
#w<-setdiff(w,sample(w,min(c(8000,length(w)))))
#obs<-obs[-w,]

#spobs<-obs[species==sp,]
#obs<-obs[species!=sp,]
#nsamp<-1000
#ran<-spobs[1:nsamp,]
#add<-st_as_sf(st_sample(region,size=nsamp,type="random"))
##add<-st_jitter(add,amount=50)
#st_geometry(add)<-"geometry"
#ran[,x:=st_coordinates(add)[,1]]
#ran[,y:=st_coordinates(add)[,2]]
#spobs[,x:=x+0+1*rnorm(nrow(spobs),0,50)]
#spobs[,y:=y+1*rnorm(nrow(spobs),0,50)]
#obs<-rbind(obs,ran,spobs)
#s<-st_as_sf(obs,coords=coords,crs=st_crs(na))
#o<-st_intersects(s,dmesh)
#o[sapply(o,function(i){length(i)==0L})] <- NA
#obs[,dmesh:=unlist(o)]

### add obs in a predictor
#varia<-"builtup_esa"
#wm<-which.max(dmeshPred[,varia])
#add<-obs[1:2,]
#xy<-st_coordinates(st_centroid(dmesh[wm,]))
#add[,x:=xy[1,1]]
#add[,y:=xy[1,2]]
#add[,dmesh:=wm]
#add[,species:=sp]
#obs<-rbind(obs,add)


### species obs
spobs<-obs[species==sp,]
#spobs<-target ### for common nighthawk

#rev(sort(table(spobs$recordedBy)))[1:5]
#plot(spobs$x[spobs$recordedBy=="Bill Chambers"],spobs$y[spobs$recordedBy=="Bill Chambers"],xlim=c(1935,1940))

nndist<-knn.dist(as.matrix(spobs[,c("x","y")]),k=1)[,1]
spobs<-spobs[nndist<=1200,]

#plot(st_geometry(st_convex_hull(st_union(dmesh))))
#plot(st_geometry(st_as_sf(spobs,coords=c("x","y"),crs=crsr)),add=TRUE)
#plot(st_geometry(region),border="red",add=TRUE)


### Drop dmesh obs info from previous run
drops<-c("nbobs","nbsp","spobs")
matches<-match(drops,names(dmesh))
if(all(!is.na(matches))){
  dmesh<-dmesh[,-matches]
}


#######################################################
### Effort info #######################################
nbobs<-obs[,.(nbobs=.N),by=dmesh]
nbsp<-obs[,.(nbsp=length(unique(species))),by=.(dmesh)]
dmesh$dmesh<-dmesh$id
dmesh<-merge(dmesh,nbobs,all.x=TRUE)
dmesh$nbobs[is.na(dmesh$nbobs)]<-0
dmesh<-merge(dmesh,nbsp,all.x=TRUE)
dmesh$nbsp[is.na(dmesh$nbsp)]<-0

#occs<-loccs[[j]]
### Species effort
temp<-spobs[,.(spobs=.N),by=.(dmesh)]
dmesh<-merge(dmesh,temp,all.x=TRUE)
dmesh$spobs[is.na(dmesh$spobs)]<-0
dmesh$pres<-as.integer(dmesh$spobs>0)
vals<-dmesh$nbobs*scales::rescale(dmesh$pres/dmesh$nbsp,to=c(1,max(dmesh$nbsp)^1))
#vals<-dmesh$nbobs*scales::rescale(dmesh$spobs/dmesh$nbobs,to=c(1,max(dmesh$nbsp)^1))
#vals<-dmesh$nbobs
vals<-ifelse(is.nan(vals) | is.infinite(vals),0,vals)
dmesh$effoccs<-vals
#k<-dmesh$spobs==0 & dmesh$nbobs>0 #############
#vals[k]<-dmesh$nbobs[k]*dmesh$nbsp[k] ###########
#dmesh$effoccs[k]<-vals[k] ############


#dmesh$nboccs<-lengths(st_intersects(dmesh,occs))
#dmesh$pres<-as.integer(dmesh$nboccs>0)
#vals<-dmesh$nbobs*scales::rescale(dmesh$pres/dmesh$nbsp,to=c(1,max(dmesh$nbsp)^1))
#vals<-ifelse(is.nan(vals) | is.infinite(vals),0,vals)
#dmesh$effoccs<-vals


#occs<-s[s$species==x,]  
#### removes points with nn > 1000 km
#nndist<-knn.dist(st_coordinates(occs),k=1)[,1]
#occs<-occs[!nndist>=1200,]
#occs

#######################################################
### Add artificial effort #############################
# adds an artificial effort value when at a certain distance bdist from concave hull 
nline<-48
effvalue<-20#50
bdist1<-500#500
bdist2<-1500#500
#buff1<-st_buffer(concaveman(st_as_sf(st_cast(st_union(st_as_sf(spobs,coords=c("x","y"),crs=st_crs(na))),"MULTIPOINT"),concavity=1)),bdist1)
#buff2<-st_buffer(concaveman(st_as_sf(st_cast(st_union(st_as_sf(spobs,coords=c("x","y"),crs=st_crs(na))),"MULTIPOINT"),concavity=1)),bdist2)
buff1<-st_union(st_buffer(st_as_sf(st_cast(st_as_sf(spobs,coords=c("x","y"),crs=st_crs(na)),"MULTIPOINT")),bdist1))
buff2<-st_union(st_buffer(st_as_sf(st_cast(st_as_sf(spobs,coords=c("x","y"),crs=st_crs(na)),"MULTIPOINT")),bdist2))
o1<-!as.logical(lengths(st_intersects(dmesh,buff1)))
o2<-!as.logical(lengths(st_intersects(dmesh,buff2)))
o<-o1 & !o2 & st_coordinates(st_centroid(st_transform(dmesh,4326)))[,2]>nline
o<-!o & o1 
dmesh$effoccs[o]<-ifelse(dmesh$effoccs[o]>0,dmesh$effoccs[o],effvalue)

#png("/data/sdm_rbq/graphics/effort_buffer.png",width=5,height=5,res=400,units="in")
#par(mar=c(0,0,0,0))
#plot(st_geometry(na),col=colmean[1],border=NA,axes=FALSE)
#naplot(lwd=0.5)
#plot(st_geometry(getobs(sp)),pch=16,cex=0.7,col=adjustcolor(colmean[95],0.5),add=TRUE)
#plot(st_geometry(dmesh[o,]),col=adjustcolor("seagreen",0.3),border=NA,add=TRUE)
#lon<-seq(-180,-50,by=1)
#lat<-rep(nline,length(lon))
#lat<-st_as_sf(data.frame(lon,lat),coords=c("lon","lat"),crs=4326) |> 
#  st_combine() |> 
#  st_cast("LINESTRING") |>
#  st_transform(st_crs(na))
#lat<-st_intersection(lat,dmesh)
#plot(st_geometry(lat),col=adjustcolor("black",0.5),lwd=3,add=TRUE)
#legend("topright",legend=c("Assumed absence",paste0("Latitude ",nline,"\u00B0N"),"Observation"),pch=c(15,NA,16),col=c(adjustcolor("seagreen",0.3),adjustcolor("black",0.95),adjustcolor(colmean[95],0.75)),pt.cex=c(1.25,1,0.9),bty="n",cex=0.9,lwd=c(NA,3,NA),seg.len=1)
#dev.off()

#dmesh$effoccs<-efftarget ### for common night hawk

########################################################
### Models #############################################

### Remove variables with very low coverage

vars<-vars_pool

if(TRUE){
  # remove interactions
  x<-unlist(lapply(grep("2|3|4|5|6",vars_pool,value=TRUE,invert=TRUE),function(k){
    x<-backScale(dmeshPred[,k],k)
    h<-hist(x[dmesh$spobs>0],breaks=seq(min(x),max(x),length.out=100),  plot=FALSE)
    #h$density
    w<-range(which(h$density>0))
    #if(sum(h$density[11:length(h$density)])==0){
    if(k=="elevation"){
      if(quantile(x[dmesh$spobs>0],0.75)>=1000){
        NULL
      }else{
        k
      }  
    }else{
      if((length(w[1]:w[2])/length(h$density))<=0.20){ # removes vars with coverage below 0.20
      #if(h$mids[w]<=0.25){ # removes vars with coverage below 0.15
        k
      }else{
        NULL
      }
    }
    # if max positive value is inferior to 5% or 10% of values, remove this variable
  }))
  if(!is.null(x)){
    x<-c(x,paste0(x,rep(2:6,each=length(x))))
    vars<-vars_pool[!vars_pool%in%x]
    vars<-unique(c("tmean",vars)) # always keep tmean as a simple if coverage is low
  }else{
    vars<-vars_pool
  }
}


### Params
f<-formula(paste("y~",paste(vars,collapse="+")))
formula<-f
sPoints<-st_as_sf(spobs,coords=c("x","y"),crs=crsr)
#explanaMesh<-explana
ppWeight<-weights
prior.range<-c(50,0.1)
prior.sigma<-c(1,0.1)
smooth<-3/2
num.threads<-1:1
blas.num.threads<-1
many<-TRUE
control.inla<-list(strategy="adaptive",int.strategy="eb",huge=TRUE) # adaptive, eb
fix<-NULL
sboffset<-"sbias"
inla.mode<-"experimental"
control.fixed<-bpriors
control.compute<-list(config=TRUE,openmp.strategy="pardiso.parallel")
verbose<-TRUE


#==============
# Basic objects
#==============
#nsPoints <- length(sPoints)
nsPoints <- nrow(sPoints)
nEdges <- Mesh$n
xy <- st_coordinates(sPoints)
colnames(xy) <- c("x","y") 

#============
# Define SPDE
#============
SPDE <- inla.spde2.pcmatern(mesh=Mesh,
                            alpha=smooth,
                            prior.range=prior.range,
                            prior.sigma=prior.sigma,
                            constr=TRUE)
  
  #====================================================
  # Rescale weights if sampling bias offset is included
  #====================================================

dmesh$weights<-dmesh$areas
  
if(!is.null(sboffset)){
    #checkpoint("Extracting sboffset")
    #e <- exact_extract(explanaMesh$X[[sboffset]], 
    #    attributes(ppWeight)$dmeshcuts, 
    #    fun = function(values, coverage){
    #    sum(values * coverage, na.rm = TRUE)
    #},progress = FALSE)
    #k <- ppWeight > 0
    #ppWeight[k] <- ppWeight[k] * ((e/ppWeight[k])/max(e/ppWeight[k]))
    #checkpoint("Done")

    k <- dmesh$areas > 0
    e <- dmesh$effoccs[k]
    dmesh$weights[k]<-dmesh$areas[k] * ((e/dmesh$areas[k])/max(e/dmesh$areas[k]))
}
  
  
#========================
# Define response objects
#========================

# Aggregate spatial data
xyDF <- as.data.frame(xy)

#spaceAgg <- mapSpecies:::aggData(sPoints, meshSpace = explanaMesh$meshSpace, meshDual = attributes(ppWeight)$dmesh)
spaceAgg <- data.frame(space = 1:nrow(dmesh), Freq = dmesh$spobs)
# Pseudo-absences are the number of edges on the mesh
# Occurences are the number of points
yPP <- spaceAgg$Freq
# weight associated to pseudo-absences (w) and occurrences (0)
ePP <- dmesh$weights[spaceAgg$space]
checkpoint("Done")

XEst<-dmeshPred[,vars,drop=FALSE]
#XEst<-apply(XEst,2,function(i){scales::rescale(i,0:1)})
XPred <- XEst


#=====================
# Fix given predictors
#=====================
if(!is.null(fix)){
    m <- match(fix, colnames(XPred))
    v <- apply(XPred[, m, drop=FALSE], 2, max, na.rm = TRUE)
    XPred[, m] <- rep(v, each = nrow(XPred))
}

#================================================
# Define projection matrix and build stack object
#================================================
# Projection matrix
ProjInfer <- inla.spde.make.A(Mesh,loc = Mesh$loc[spaceAgg$space,])
IDSpace <- inla.spde.make.index('i', nEdges)

# Build stack objects
StackEst <- inla.stack(data = list(y = yPP, e = ePP), 
                            A = list(1, ProjInfer), 
                            effects = list(c(list(Intercept = 1), 
                                            asplit(XEst, 2)),
                                        IDSpace), 
                            tag = "est")
    
StackPred <- inla.stack(data = list(y = NA, e = NA),
                            A = list(1, ProjInfer), 
                            effects = list(c(list(Intercept = 1), 
                                            asplit(XPred, 2)), 
                                            IDSpace),
                            tag = "pred")

Stack <- inla.stack(StackEst, StackPred)
  

X <- paste(colnames(XEst),collapse=" + ")
fixed <- paste("y ~ 0 + Intercept +",X)
formule <- formula(paste(fixed,"f(i, model=SPDE)", sep=" + "))

if(TRUE){
  # build constraints
  XX = cbind(rep(1, nrow(XEst)), XEst)
  Q = qr.Q(qr(XX))
  AA <- as.matrix(t(Q)%*%ProjInfer)
  ee <- rep(0, ncol(XX))
  formule <- formula(paste(fixed,
                                "f(i, model=SPDE, extraconstr = list(A = AA, e = ee))",
                                sep=" + "))
}




verbose<-TRUE  
checkpoint("Running model")
model <- inla(formule, 
              family = "poisson", 
              data = inla.stack.data(Stack),
              control.predictor = list(A = inla.stack.A(Stack),link = 1),
              E = inla.stack.data(Stack)$e, 
              num.threads=4:4,
              blas.num.threads=4,
              control.inla=list(strategy="adaptive",int.strategy="eb",huge=TRUE,
              control.vb=list(enable=TRUE, verbose=verbose)),
              inla.mode="experimental",
              control.fixed=bpriors,
              control.compute=list(config=TRUE,openmp.strategy="pardiso"),
              verbose=verbose
             )


nameRes <- names(model)
  
#prog(message = sprintf("Model done %s", species[j]))

#=============
# Return model
#=============

attributes(model) <- list(formula = formula,
                            sPoints = sPoints,
                            XEst = XEst,
                            XPred = XPred,
                            meshSpace = Mesh,
                            Stack = Stack)

names(model) <- nameRes

checkpoint("Plotting marginal effects")
plot_marginal_effects(m=model,sp=sp,dmesh=dmesh)

class(model) <- "ppSpace"

#mpp<-model

##################################################################
### Ouputs dependent on models ###################################

gb<-colorRampPalette(c("grey94","lightsteelblue1"))(100)[5]

colo<-list(
  #mean=colorRampPalette(c("grey90","steelblue","steelblue2","gold2","tomato2","red4"))(200),
  #mean=colorRampPalette(c("grey90","tomato","darkred","black"))(200)[1:200],
  mean=colorRampPalette(c(gb,"orangered3","red4"))(200)[1:200],
  #mean=c("snow2",unname(palette.colors()[7]),"darkred","darkred","grey10"),
  sd=rev(magma(200))[1:175]
)

colmean<-c("snow2",unname(palette.colors()[7]),"darkred","darkred","grey10")

colmean<-colo.scale(1:200,c("#EEEFF0","#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00","darkorange")) #https://observablehq.com/@d3/radial-stacked-bar-chart


checkpoint("Chull reach asymptote?")
rh<-rangeHull(sPoints,species=sp,breaks=200)
checkpoint("Plot stabilization")
stabHull(sp,rh=rh)
checkpoint("Plot done")
rh<-rh$reach


checkpoint("Mapping space for pixelate and animate")
sdm<-mapSpace(modelSpace=model,dims=round(1*c(1,1)*dims/4,0),sPoly=NULL,sample=TRUE)#[[type]]
sdm<-exp(sdm[[grep("sample",names(sdm))]])
sdm<-crop(sdm,vect(na))
sdm<-mask(sdm,vect(na))

plot(1,1);coloScale(vals=1:200)
checkpoint("Pixelate")
pixelate(sdm=sdm,sp=sp)
checkpoint("Animate")
animate(sdm=sdm,sp=sp)


checkpoint("Mapping distribution")
png(paste0("/data/sdm_rbq/sdms/","birds_",gsub(" ","_",sp),".png"),units="in",width=10,height=8,res=200)
m<-list(mpp=model)
par(mfrow=n2mfrow(length(m)),mar=c(0,0,0,0))
for(i in 1:length(m)){
  type<-"mean"
  cols<-colo[[type]]
  sdm<-mapSpace(modelSpace=m[[i]],dims=round(1*c(1,1)*dims,0),sPoly=NULL,sample=FALSE)#[[type]]
  sdm <- mask(sdm,vect(na)) #spacePoly
  #buff<-st_as_sf(st_buffer(st_convex_hull(st_union(occs)),dist=1000))
  #sdm <- mask(sdm,vect(buff)) #spacePoly
  #sdm <- mask(sdm,vect(nalakes),inverse=TRUE) #spacePoly
  #names(sdm)<-paste(sp,names(m)[i])
  plot(st_geometry(na),border=NA,col="white")
  plot(sdm[[type]],add=FALSE, col = cols, axes = FALSE, main = "")
  #plot(st_geometry(em[em$species==species[j],]),col=gray(0,0.2),border=NA,add=TRUE)
  plot(st_geometry(na),add=TRUE,lwd=0.75,border=gray(1,0.99))
  plot(st_geometry(nalakes),col="white",border=adjustcolor("black",0.15),lwd=0.5,add=TRUE)
  plot(st_geometry(st_union(na)),lwd=0.5,border=gray(0,0.45),add=TRUE)
  plot(st_geometry(sPoints),pch=16,col=adjustcolor("black",0.15),cex=0.5,add=TRUE)
  ### contour
  ID<-inla.stack.index(attributes(m[[i]])$Stack,tag="est")$data
  if(type=="mean"){
    #lim<-rasterToContour(raster(sdm[[type]]),levels=max(values(sdm[[type]]),na.rm=TRUE)*0.20)
    #plot(lim,add=TRUE,col=alpha("black",0.65),lwd=1)
  }
  crs(sdm)<-crsr
  period<-format(as.Date(paste("2000",be,sep="-")),format="%b-%d")
  lab<-paste(paste(sp,""),paste(period,collapse=" / "),paste("n =",nrow(spobs)),paste("reach =",round(rh,2)),sep="   ")
  mtext(side=3,line=-1,adj=0.05,text=lab,font=2,cex=1.4,col=gray(0,0.3))
  #plot(attributes(weight)$dmesh,add=TRUE,border=gray(1,0.25),lwd=0.2)
  #mf_scale()
}
par(mfrow=c(1,1))
dev.off()


writeRaster(sdm, filename=paste0("/data/sdm_rbq/rasters/",gsub(" ","_",sp),"_birds.tif"), overwrite=TRUE)
 
time<-Sys.time()
attr(time,"tzone")<-"Indian/Reunion"
res<-list(species=sp,mpp=model,spobs=spobs,dmesh=dmesh,n=nrow(spobs),reach=round(rh,5),date=time,predictors=paste(vars,collapse="+"),range=model$summary.hyperpar[1,1],sd=model$summary.hyperpar[2,1])

### write results
results<-data.frame(
    species=res$species,
    date=res$date,
    n=res$n,
    reach=res$reach,
    predictors=res$predictors,
    range=res$range,
    sd=res$sd,
    pearson=NA,ratio
    spearman=NA,
    I=NA,
    D=NA,
    hullarea=st_union(sPoints) |> st_convex_hull() |> st_area() |> as.numeric(),
    family=d$family[match(res$species,d$species)],
    familyname=d$familyname[match(res$species,d$species)],
    order=d$order[match(res$species,d$species)],
    ordername=d$ordername[match(res$species,d$species)],
    max=NA,
    hullratio=NA

  )
write.table(results,file="/data/sdm_rbq/graphics/mapSpeciesres.csv",row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)

print(j) 
res

#prog(message = sprintf("Results done %s", species[j]))

}, future.packages = "data.table", future.seed=TRUE)
plan(sequential)


### clear old results
res<-rev(sort(list.files("/data/sdm_rbq/graphics",pattern="mapSpeciesres",full=TRUE)))[1]
df<-read.table(res,sep=",",header=TRUE)
df<-df[rev(order(df$date)),]
df<-df[!duplicated(df$species),]
#tz<-as.POSIXct(df$date)
#attr(tz,"tzone")<-"Indian/Reunion"
#df$date<-tz
write.table(df,file="/data/sdm_rbq/graphics/mapSpeciesres.csv",row.names=FALSE,sep=",",append=FALSE)



#df<-df[,c("species", "date", "n", "reach", "predictors", "range", "sd", "pearson","spearman", "I", "D")]
#lsps<-gsub("_"," ",gsub("_birds.tif","",basename(list.files("/data/sdm_rbq/rasters",full=TRUE)[which(file.info(list.files("/data/sdm_rbq/rasters",full=TRUE))$mtime>as.POSIXct("2022-11-16"))])))
#lsps<-lsps[!lsps%in%df$species]

#results<-data.frame(
#    species=lsps,
#    date=Sys.time(),
#    n=NA,
#    reach=NA,
#    predictors=NA,
#    pearson=NA,
#    spearman=NA,
#    I=NA,
#    D=NA
#  )
#write.table(results,file="/data/sdm_rbq/graphics/mapSpeciesres.csv",row.names=FALSE,sep=",",append=FALSE)


species
#i<-1
#sp<-species[i]
#mpp<-ml[[i]]$mpp
#sdm<-rast(paste0("/data/sdm_rbq/rasters/",gsub(" ","_",species[i]),"_birds.tif"))
#occs<-st_as_sf(ml[[i]]$spobs,coords=c("x","y"),crs=crsr)
#dmesh<-ml[[i]]$dmesh
#sdmf<-mapSpace(modelSpace=mpp,dims=round(1*c(1,1)*dims/4,0),sPoly=NULL,sample=TRUE)#[[type]]
#system(paste("code",paste0("/data/sdm_rbq/sdms/",paste0("birds_",gsub(" ","_",species[i]),".png"))), wait=FALSE)

head(df,10)

#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################

if(FALSE){


library(spatstat.core)
data(redwood)
u<-lgcp.estpcf(redwood,c(var=1, scale=0.1))
u
plot(u)



#rh<-rangeHull(occs,species=occs$species[1],breaks=200)$reach
#rh

dmesh$pred<-identity(mpp$summary.fitted.values[1:nrow(dmesh),"mean"])
dmesh$pred<-mpp$summary.random$i[,"mean"]
#plot(dmesh["pred"],border=NA,nbreaks=200,logz=FALSE,reset=FALSE)
#plot(st_geometry(buff),add=TRUE);plot(st_geometry(na),add=TRUE)

#plot(dmesh["spobs"],border=NA,nbreaks=200,logz=TRUE,reset=FALSE)

res<-as.data.frame(mpp$summary.fitted.values[1:nrow(dmesh),])
res<-as.data.frame(mpp$summary.random$i)
res<-cbind(res,st_drop_geometry(dmesh[,c("areas","nbobs","spobs","pres","effoccs")]))
w<-which.max(res$sd)
#w<-identity(order(res$mean))[1:1000]
res[w,]



w<-which.max(dmesh$pred)

tampon<-st_buffer(dmesh[w,],dist=500)
reg<-st_intersection(dmesh,tampon)
plot(reg["pred"],reset=FALSE,nbreaks=200)
plot(st_geometry(tampon),add=TRUE)
plot(st_geometry(dmesh[w,]),add=TRUE)
plot(st_geometry(na),add=TRUE)
plot(st_geometry(occs),add=TRUE,cex=0.5,col=adjustcolor("black",0.5))
#plot(dmesh["pred"],key.pos=NULL,add=TRUE)


library(glmmTMB)
library(visreg)
dmesh$y<-dmesh$spobs
dat<-cbind(as.data.frame(dmeshPred),st_drop_geometry(dmesh))
dat<-dat[dat$nbobs>0,]
fit<-glmmTMB(y ~ poly(tmean,2)+offset(log(effoccs)),data=dat,family="poisson")
par(mfrow=n2mfrow(1))
v<-visreg(fit,scale="response",partial=TRUE)
dmesh$pred<-predict(fit,newdata=data.frame(dmeshPred,effoccs=100),type="response")
plot(dmesh["pred"],border=NA,nbreaks=200,logz=FALSE,reset=FALSE,pal=magma)
plot(st_geometry(na),add=TRUE)



plot(dmesh["effoccs"],border=NA,nbreaks=200,logz=TRUE,reset=FALSE)
plot(st_geometry(na),add=TRUE)
plot(st_geometry(dmesh[dmesh$dmesh==980,]),add=TRUE,col="red")




res1<-as.data.frame(mpp$summary.fitted.values[1:nrow(dmesh),])

res2<-as.data.frame(mpp$summary.linear.predictor[1:nrow(dmesh),])

c(res1[10,"mean"],exp(res2[10,"mean"]))


dmesh$vals<-identity(dmesh$effoccs)
#dmesh$vals<-dmesh$nbsp
plot(dmesh["vals"],border=NA,reset=FALSE)
plot(st_geometry(occs),cex=0.1,add=TRUE)


dmesh$weights<-dmesh$areas
k <- dmesh$areas > 0
e <- dmesh$effoccs[k]
dmesh$weights[k]<-dmesh$areas[k] * ((e/dmesh$areas[k])/max(e/dmesh$areas[k]))

range(dmesh$weights)


gb<-colorRampPalette(c("grey94","lightsteelblue1"))(100)[5]
par(mar=c(0,0,0,0))
plot(st_geometry(na),col=gb,border="white",lwd=0.75)


plot(dmesh["effoccs"],border=NA,logz=TRUE)

plot(dmesh["effoccs"][dmesh$spobs==0,],border=NA,logz=TRUE,reset=FALSE)
plot(st_geometry(na),add=TRUE)





mapSpace2<-function (modelSpace, locs = NULL, dims, sPoly = NULL, sample = FALSE, nsamples = 100) 
{
    valsPred <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", 
        "mode")
    valsLink <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", 
        "mode")
    valsSpat <- c("mean", "sd")
    if (sample) {
        valsSamp <- paste0("sample", formatC(1:nsamples, width = nchar(nsamples), 
            flag = 0))
    }
    else {
        valsSamp <- NULL
    }
    if (is.null(sPoly)) {
        mapBasis <- inla.mesh.projector(attributes(modelSpace)$mesh, loc = locs, dims = dims, crs = attributes(modelSpace)$mesh$crs$crs)
    }
    else {
        mapBasis <- inla.mesh.projector(attributes(modelSpace)$mesh, loc = locs, dims = dims, xlim = c(xmin(sPoly), xmax(sPoly)), 
            ylim = c(ymin(sPoly), ymax(sPoly)), crs = attributes(modelSpace)$mesh$crs$crs)
    }
    mapSpat <- as.matrix(inla.mesh.project(mapBasis, modelSpace$summary.random[["i"]][, 
        valsSpat]))
    ID <- inla.stack.index(attributes(modelSpace)$Stack, tag = "pred")$data
    mapPred <- as.matrix(inla.mesh.project(mapBasis, modelSpace$summary.fitted.values[ID, 
        valsPred]))
    mapLink <- as.matrix(inla.mesh.project(mapBasis, modelSpace$summary.linear.predictor[ID, 
        valsLink]))
    if (sample) {
        class(modelSpace) <- "inla"
        samps <- inla.posterior.sample(nsamples, modelSpace)
        samps <- lapply(samps, function(i) {
            i$latent
        })
        samps <- do.call("cbind", samps)
        vals <- samps[grep("i:", row.names(samps)), ]
        vals <- samps[1:nrow(vals), ]
        mapSamp <- as.matrix(inla.mesh.project(mapBasis, vals))
        mat <- cbind(mapPred, mapLink, mapSpat, mapSamp)
    }
    else {
        mat <- cbind(mapPred, mapLink, mapSpat)
    }
    a <- array(rev(as.vector(mat)), dim = c(dims[1], dims[2], 
        ncol(mat)))
    a <- apply(a, 3, t, simplify = FALSE)
    a <- lapply(a, function(i) {
        i[, ncol(i):1]
    })
    a <- simplify2array(a)
    mapRaster <- rast(a[, , dim(a)[3]:1])
    x<-if(!is.null(locs)){locs[,1]}else{mapBasis$x}
    y<-if(!is.null(locs)){locs[,2]}else{mapBasis$y}
    ext(mapRaster) <- c(xmin = min(x), xmax = max(x), 
        ymin = min(y), ymax = max(y))
    names(mapRaster) <- c(valsPred, paste0("link", valsLink), 
        paste0("space", valsSpat), valsSamp)
    mapRaster
}

xy<-aggregate(r[[1]],10)
xy<-xyFromCell(xy,1:ncell(xy))
colnames(xy)<-c("x","y")
sdmf<-mapSpace2(modelSpace=mpp,locs=NULL,dims=c(100,100),sPoly=NULL,sample=FALSE)#[[type]]



        mapBasis <- inla.mesh.projector(Mesh, loc = NULL, dims = dims, crs = attributes(mpp)$mesh$crs$crs)






plot(sdm$mean,col=FRutils::colo.scale(1:200,colmean))

dmesh$pred<-identity(mpp$summary.fitted.values[1:nrow(dmesh),"mean"])



i<-rev(order(dmesh$pred))[1:10]

lcc<-"mixed"
dmesh$predictor<-dmeshPred[,lcc]
hist(backTransform(dmesh$predictor[dmesh$spobs>0],lcc),breaks=seq(0,1,by=0.025),border=FALSE)

# plot(dmesh["pred"],pal=magma,nbreaks=200,border=NA,reset=FALSE)
plot(dmesh["predictor"],pal=magma,nbreaks=200,border=NA,reset=FALSE)
plot(st_geometry(dmesh["predictor"][i,]),border="red",add=TRUE)
plot(st_geometry(na),border=adjustcolor("white",0.5),add=TRUE)

#plot(dmesh["pred"],pal=magma,nbreaks=200,border=NA,add=TRUE)
plot(dmesh["predictor"],pal=magma,nbreaks=200,border=NA,reset=FALSE)
plot(st_geometry(na),border=adjustcolor("white",0.5),add=TRUE)

}




if(FALSE){

remotes::install_github("frousseu/ewlgcpSDM",auth_token=readLines("/home/rouf1703/github_token_pose"),upgrade="never")

flex_buffer<-function(
    obs,
    nline=48,
    dist=c(500,1500)
){
  # https://r-spatial.github.io/sf/articles/sf7.html#buffers-1

  radius<-s2_earth_radius_meters()*pi*2*0.25*((90-nline)/90)
  npole<-s2_buffer_cells(as_s2_geography("POINT(-100 90)"),distance=radius,max_cells=10000) # visible half
  g<-as_s2_geography(TRUE)
  north<-s2_intersection(npole, g)
  north<-st_transform(st_as_sfc(north),st_crs(na))



  buff1<-st_union(st_buffer(obs[lengths(st_intersects(obs,north))==0L,],dist[1]))
  buff2<-st_union(st_buffer(obs[lengths(st_intersects(obs,north))>0L,],dist[2]))

  bnorth<-st_cast(st_intersection(buff2,north),"POLYGON")
  o<-st_intersects(bnorth,na)
  bnorth<-bnorth[lengths(o)>0]
  #bsouth<-st_cast(st_difference(buff1,north),"POLYGON")
  bsouth<-st_cast(buff1,"POLYGON")
  o<-st_intersects(bsouth,na)
  bsouth<-bsouth[lengths(o)>0]
  buff<-st_union(bnorth,bsouth) |> st_union()
  #plot(st_geometry(na))
  #plot(st_geometry(buff),col=adjustcolor("black",0.15),border=NA,add=TRUE)
  #plot(st_geometry(obs),add=TRUE)
  #plot(st_geometry(bnorth),col="blue",add=TRUE)
  #plot(st_geometry(bsouth),col="red",add=TRUE)
  #lpols<-lapply(buff, function(x) x[1])
  #buff<-st_multipolygon(lpols) # remove holes https://github.com/r-spatial/sf/issues/609
  buff
}


library(ewlgcpSDM)
library(predicts)
library(s2)

sp<-"Spizella pusilla" 

be<-unlist(ed[match(sp,ed$scientific_name),c("start","end")])
### all obs for period
if(be[1]<be[2]){
  background<-d[md>=be[1] & md<=be[2],]
}else{
  # Tyrannus savanna breeds from nov to jan
  dr<-substr(seq.Date(as.Date(paste0("2007-",be[1])),as.Date(paste0("2008-",be[2])),by=1),6,10)
  background<-d[md%in%dr,]
}
background<-unique(background,by=c("recordedBy","species","cell"))
obs<-background[species==sp,]

obs<-st_as_sf(obs,coords=c("x","y"),crs=st_crs(na))
background<-st_as_sf(background,coords=c("x","y"),crs=st_crs(na))
#buff<-st_buffer(st_combine(obs),500)
#buff<-st_union(st_buffer(obs,500))
buff<-flex_buffer(obs,nline=48,dist=c(500,1500))


tb1<-background[sample(1:nrow(background),100000),]
tb2<-st_as_sf(st_sample(st_union(st_difference(na,buff)),50000))
st_geometry(tb2)<-"geometry"
tb<-rbind(tb1[,names(tb2)],tb2)
#plot(st_geometry(tb),add=TRUE)


# fit model
#mep<-aggregate(predictors[[c("tmean","trange","prec","elevation","truggedness","latitude","longitude",grep("_esa",names(predictors),value=TRUE))]],8,na.rm=TRUE)
me <- MaxEnt(mep,vect(obs),vect(tb),args=c("-q","-l","-p"))
#me <- MaxEnt(mep,vect(obs),args=c("-p","-h"))
mer <- predict(me, mep)
mer <- mask(mer,vect(na))
plot(mer,mar=c(0,0,0,8),plg=list(size=c(0.8,1)),col=colmean)
naplot()
#plot(st_geometry(obs),pch=1,cex=1,col=adjustcolor("black",0.5),add=TRUE)

i<-gsub(" ","_",sp)
r1<-raster(file.path("/data/sdm_rbq/rasters",paste0(i,"_birds.tif")),layer=1)
r2<-raster(file.path("/data/predictors_sdm/expert_maps/eBird/abundance",paste0(i,"_ebird2.tif")),layer=1)
r3<-resample(raster(mer),r2)
nicheOverlap(r1,r2,stat="I")
nicheOverlap(r3,r2,stat="I")
#nicheOverlap(r1,r3,stat="I")
#r4<-mean(c(rast(r1)/global(rast(r1),"max",na.rm=TRUE)[1,1],rast(r3)/global(rast(r3),"max",na.rm=TRUE)[1,1]));plot(r4,col=colmean);#nicheOverlap(raster(r4),r2,stat="I")
#r4<-rast(r1)/global(rast(r1),"max",na.rm=TRUE)[1,1]*rast(r3)/global(rast(r3),"max",na.rm=TRUE)[1,1];plot(r4,col=colmean);#nicheOverlap(raster(r4),r2,stat="I")

r11<-rast(r1)
r11<-r11/global(r11,"max",na.rm=TRUE)[1,1]
th<-0.1
r11[r11<th]<-NA
r11[r11>=th]<-1
v <- as.polygons(r11)
v<-st_as_sf(v)

r11<-rast(r1)
r11<-r11/global(r11,"max",na.rm=TRUE)[1,1]
r11<-mask(r11,vect(v),inverse=TRUE,updatevalue=1)
r4<-r11*rast(r3)

r11
r4<-rast(r1)/global(rast(r1),"max",na.rm=TRUE)[1,1]*rast(r3)/global(rast(r3),"max",na.rm=TRUE)[1,1]


plot(rast(r1),col=colmean)
naplot()
plot(st_geometry(v),col=adjustcolor("black",0.15),add=TRUE)

r33<-rast(r3)
r33<-mask(r33,vect(v),inverse=FALSE,updatevalue=0)
r33<-mask(r33,vect(na))
plot(r33,col=colmean)
naplot()
plot(st_geometry(v),col=adjustcolor("black",0.15),add=TRUE)
#nicheOverlap(raster(r33),r2,stat="I")




plan(multicore,workers=5)
params<-dmesh_mesh(Mesh)
params<-dmesh_weights(params,region)
params<-dmesh_predictors(params,r)
params<-dmesh_effort(params,obs=obs,background=background,buffer=buff,adjust=TRUE)
plan(sequential)


m<-ewlgcp(
  formula=formula(paste("y~",paste0(setdiff(vars_pool,"logdistance"),collapse="+"))),
  dmesh=params,
  effort = TRUE,
  adjust = FALSE,
  buffer = TRUE,
  orthogonal = TRUE,
  prior.beta = NULL,
  prior.range = c(50,0.1),
  prior.sigma = c(1,0.1),
  smooth = 3/2
)

sdm<-ewlgcpSDM::map(model=m,
         dmesh=params,
         dims=c(1500,1500),
         region = region,
         sample = FALSE,
         nsamples = 100
     )

sdm<-mask(sdm,vect(na))


plot(sdm$mean,col=colmean,plg=list(size=c(0.8,1)))
naplot()
plot(st_geometry(obs),add=TRUE)



#r5<-raster(resample(sdm$mean,rast(r1)))


}


dm<-params$dmesh

o<-st_intersects(background,dm)
splist<-data.frame(species=background$species,id=dm$id[unlist(o)])
res<-aggregate(species~id,
               data=splist,FUN=function(i){
                     length(unique(i))
               }
)
nsp<-res$species[match(dm$id,res$id)]
nsp[is.na(nsp)]<-0
pres<-as.integer(nobs>0)
vals<-nobs*rescale_ab(pres/nsp,a=1,b=max(nsp))
vals<-ifelse(is.nan(vals) | is.infinite(vals),0,vals)




mepp<-wrap(mep)


cl<-makeCluster(20)
registerDoParallel(cl)
lsp<-sample(df$species[which(!is.na(df$I) & df$reach>=0.85)],80)
comp<-foreach(sp=lsp,.packages=c("dismo","raster","terra","predicts","sf","data.table","s2")) %dopar% {

be<-unlist(ed[match(sp,ed$scientific_name),c("start","end")])
### all obs for period
if(be[1]<be[2]){
  background<-d[md>=be[1] & md<=be[2],]
}else{
  # Tyrannus savanna breeds from nov to jan
  dr<-substr(seq.Date(as.Date(paste0("2007-",be[1])),as.Date(paste0("2008-",be[2])),by=1),6,10)
  background<-d[md%in%dr,]
}
background<-unique(background,by=c("recordedBy","species","cell"))
obs<-background[species==sp,]

obs<-st_as_sf(obs,coords=c("x","y"),crs=st_crs(na))
background<-st_as_sf(background,coords=c("x","y"),crs=st_crs(na))
#buff<-st_buffer(st_combine(obs),500)
#buff<-st_union(st_buffer(obs,500))
buff<-flex_buffer(obs,nline=48,dist=c(500,1500))


tb1<-background[sample(1:nrow(background),min(100000,nrow(background))),]
tb2<-st_as_sf(st_sample(st_union(st_difference(na,buff)),50000))
st_geometry(tb2)<-"geometry"
tb<-rbind(tb1[,names(tb2)],tb2)
#plot(st_geometry(tb),add=TRUE)


# fit model
#mep<-aggregate(predictors[[c("tmean","trange","prec","elevation","truggedness","latitude","longitude",grep("_esa",names(predictors),value=TRUE))]],8,na.rm=TRUE)
me <- MaxEnt(rast(mepp),vect(obs),vect(tb),args=c("-q","-l","-p"))
#me <- MaxEnt(mep,vect(obs),args=c("-p","-h"))
mer <- predict(me, rast(mepp))
#mer <- mask(mer,vect(na))
#plot(mer,mar=c(0,0,0,8),plg=list(size=c(0.8,1)),col=colmean)
#naplot()
#plot(st_geometry(obs),pch=1,cex=1,col=adjustcolor("black",0.5),add=TRUE)

i<-gsub(" ","_",sp)
r1<-raster(file.path("/data/sdm_rbq/rasters",paste0(i,"_birds.tif")),layer=1)
r2<-raster(file.path("/data/predictors_sdm/expert_maps/eBird/abundance",paste0(i,"_ebird2.tif")),layer=1)
r3<-resample(raster(mer),r2)
comp<-c(
  nicheOverlap(r1,r2,stat="I"),
  nicheOverlap(r3,r2,stat="I")
)
}


stopCluster(cl)
mmcomp<-do.call("rbind",comp)
mmcomp
colMeans(mmcomp)