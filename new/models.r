 
if(FALSE){
  source("/data/sdm_rbq/functions.r")
  source("/data/sdm_rbq/data.r")
  source("/data/sdm_rbq/parameters.r")
  source("/data/sdm_rbq/inputs.r")
}

library(terra)
library(sf)
library(concaveman)
library(future)
library(future.apply)
library(future.callr)
library(viridisLite)


### This runs models for each species and produces necessary model outputs
# The number of cores used is adjusted according to the number of species

inla.setOption(pardiso.license = "/home/rouf1703/pardiso.lic")

dims<-dim(r)[1:2]
crsr<-crs(r)

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

plan(multisession,workers=min(c(3,length(species))))
cores<-nbrOfWorkers() # get nbr of workers from the chosen plan
options(future.globals.maxSize = 5000 * 1024 ^ 2)


ml<-future_lapply(seq_along(species),function(j){

#sink("/data/sdm_rbq/graphics/output.txt", append=TRUE)
sp<-species[j] 
be<-unlist(ed[match(sp,ed$scientific_name),c("start","end")])

### all obs for period
obs<-d[md>=be[1] & md<=be[2],]

### remove duplicate obs
obs<-unique(obs,by=c("recordedBy","species","cell"))

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


### species obs
spobs<-obs[species==sp,]

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
effvalue<-50
bdist<-500 #500
buff<-st_buffer(concaveman(st_as_sf(st_cast(st_union(st_as_sf(spobs,coords=c("x","y"),crs=st_crs(na))),"MULTIPOINT"),concavity=2)),bdist)
o<-!as.logical(lengths(st_intersects(dmesh,buff)))
dmesh$effoccs[o]<-ifelse(dmesh$effoccs[o]>0,dmesh$effoccs[o],effvalue)

########################################################
### Models #############################################

### Remove variables with very low coverage

if(TRUE){
  x<-unlist(lapply(grep("2|tmean",vars_pool,value=TRUE,invert=TRUE),function(k){
    x<-backTransform(dmeshPred[,k],k)
    h<-hist(x[dmesh$spobs>0],breaks=seq(min(x),max(x),length.out=100),  plot=FALSE)
    h$density
    if(sum(h$density[8:length(h$density)])==0){
      k
    }else{
      NULL
    }
    # if max positive value is inferior to 5% or 10% of values, remove this variable
  }))
  if(!is.null(x)){
    vars<-vars_pool[!vars_pool%in%x]
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
prior.range<-c(50,0.01)
prior.sigma<-c(1,0.01)
smooth<-3/2
num.threads<-4:4
blas.num.threads<-4
many<-TRUE
control.inla<-list(strategy="adaptive",int.strategy="eb",huge=TRUE) # adaptive, eb
fix<-NULL
sboffset<-"sbias"
inla.mode<-"experimental"
control.fixed<-bpriors
control.compute<-list(config=TRUE,openmp.strategy="pardiso.parallel")
verbose<-TRUE


# Function that prints the different steps
checkpoint<-function(msg=""){
  if(verbose){
    cat(paste(Sys.time(),msg,sep=" - "),"\n")
  }
}



  
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




  
checkpoint("Running model")
model <- inla(formule, 
              family = "poisson", 
              data = inla.stack.data(Stack),
              control.predictor = list(A = inla.stack.A(Stack),link = 1),
              E = inla.stack.data(Stack)$e, 
              num.threads=4:4,
              blas.num.threads=4,
              control.inla=list(strategy="adaptive",int.strategy="eb",huge=TRUE,
              control.vb=list(enable=TRUE, verbose=TRUE)),
              inla.mode="experimental",
              control.fixed=bpriors,
              control.compute=list(config=TRUE,openmp.strategy="pardiso.parallel"),
              verbose=TRUE
             )


nameRes <- names(model)
  
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
class(model) <- "ppSpace"

mpp<-model

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

checkpoint("Chull reach asymptote?")

rh<-rangeHull(sPoints,species=sp,breaks=200)$reach

checkpoint("Mapping distribution")

png(paste0("/data/sdm_rbq/plots/","birds_",gsub(" ","_",sp),".png"),units="in",width=10,height=8,res=200)
m<-list(mpp=mpp)
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
  plot(sdm[[type]],add=FALSE, col = cols, axes = FALSE, box = FALSE, main = "")
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
  writeRaster(sdm, filename=paste0("/data/sdm_rbq/rasters/",gsub(" ","_",sp),"_birds.tif"), overwrite=TRUE)
}
par(mfrow=c(1,1))
dev.off()
 
print(j) 

list(mpp=mpp,spobs=spobs,dmesh=dmesh)
}, future.packages = "data.table")
plan(sequential)


species
i<-1
mpp<-ml[[i]]$mpp
sdm<-rast(paste0("/data/sdm_rbq/rasters/",gsub(" ","_",species[i]),"_birds.tif"))
occs<-st_as_sf(ml[[i]]$spobs,coords=c("x","y"),crs=crsr)
dmesh<-ml[[i]]$dmesh
sdmf<-mapSpace(modelSpace=mpp,dims=round(1*c(1,1)*dims/4,0),sPoly=NULL,sample=FALSE)#[[type]]


plot(mask(mean(exp(sdmf[[grep("sample",names(sdmf))]])),vect(na)))


#####################################
#####################################
#####################################
#####################################
#####################################

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

