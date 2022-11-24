
library(mapSpecies)
library(sf)
library(spatstat)
library(fields)
library(terra)
library(scales)
library(raster)
library(terra)
library(RandomFields)

grDevices:::windows.options(record=TRUE)

inla.mesh2sp <- function(mesh) {
  crs <- inla.CRS(inla.CRSargs(mesh$crs))
  isgeocentric <- identical(inla.as.list.CRS(crs)[["proj"]], "geocent")
  if (isgeocentric || (mesh$manifold == "S2")) {
    stop(paste0(
      "'sp' doesn't support storing polygons in geocentric coordinates.\n",
      "Convert to a map projection with inla.spTransform() before
calling inla.mesh2sp()."))
  }
  
  triangles <- SpatialPolygonsDataFrame(
    Sr = SpatialPolygons(lapply(
      1:nrow(mesh$graph$tv),
      function(x) {
        tv <- mesh$graph$tv[x, , drop = TRUE]
        Polygons(list(Polygon(mesh$loc[tv[c(1, 3, 2, 1)],
                                       1:2,
                                       drop = FALSE])),
                 ID = x)
      }
    ),
    proj4string = crs
    ),
    data = as.data.frame(mesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
    match.ID = FALSE
  )
  vertices <- SpatialPoints(mesh$loc[, 1:2, drop = FALSE], proj4string = crs)
  
  list(triangles = triangles, vertices = vertices)
}



###############################################
### rLGCP
###############################################

cat("\014")
seed<-sample(1:1000,1)
#seed<-285
set.seed(seed)
print(seed)

colo<-colorRampPalette(c("grey85","steelblue4","steelblue2","gold2","tomato2","red4"))(200)
#colo<-turbo(200)
prj<-"+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

inla.setOption(inla.mode="experimental")

###############################################
### Simulation parameters
###############################################



seed_intensity<-TRUE
seed_effort<-TRUE
seed_thinning<-TRUE  

win<-c(0,1000,0,1000)
o<-owin(win[1:2],win[3:4])
params<-list(
  
  ### Simulated LGCP
  nu=c(2), # matern smoothness parameters 
  rangex=c(500), # range of spatial field
  #kappa=c(sqrt(8*nu)/rangex), # matern smoothness parameter (2) and range of spatial field 
  b0=c(-5,-5,-5,-5,-5), # intercept (determines the nb of points generated
  b1=c(0.002), # beta coefficient on X gradient
  b2=c(-0.000), # beta coefficient on Y gradient
  #b3=c(-0.00001), # beta coefficient on X gradient
  #b4=c(-0.00001), # beta coefficient on Y gradient
  sigma=sqrt(c(0.000005,0.000005,0.5,0.5,0.5)), # sd of spatial field
  
  ### Simulated effort
  evar=c(0.01,0.5,1,2,3), # variance of effort spatial field
  erange=c(30), # range of effort spatial field
  emean=c(2,1,-1,-2,-3), # mean of effort spatial field
  eb1=c(-1), # strength of gradient in x
  eb2=c(-1), # strength of gradient in y
  eexp=c(-0.01,-0.01,-0.5,-1,-2) # exponential decrease in effort gradient
)

params<-lapply(params,function(i){rep(i,length.out=max(sapply(params,length)))})
params<-lapply(params,"[",c(2))

#layout(matrix(seq_len(5*length(params[[1]])),nrow=5,byrow=FALSE))

#png("plots.png",width=12,height=10,units="in",res=200)

png("C:/Users/God/Downloads/simulated_effort.png",width=19,height=length(params[[1]])*(10/4),units="in",res=300)

par(mfcol=c(length(params[[1]]),6),oma=c(0,2,0.25,3),mar=c(0.5,0.5,1,0.5))

ml<-vector(mode="list",length=length(params[[1]]))
for(k in seq_along(params[[1]])){
  
  if(seed_intensity){set.seed(seed)}
  pp<-rLGCP("matern",mu=function(x,y){params$b0[k]+params$b1[k]*x+params$b2[k]*y},var=params$sigma[k]^2,scale=1/c(sqrt(8  *params$nu[k])/params$rangex[k]),nu=params$nu[k],win=o,dimyx=c(300,300))
  #pp<-rLGCP("matern",mu=function(x,y){params$b0[k]+params$b1[k]*x+params$b3[k]*x^2+params$b2[k]*y+params$b4[k]*y^2},var=params$sigma[k]^2,scale=1/c(sqrt(8  *params$nu[k])/params$rangex[k]),nu=params$nu[k],win=o,dimyx=c(300,300))
  cutoff<-500000
  if(pp$n>cutoff){
    pn<-pp$n
    rm(pp)
    stop(paste("Too many points generated",pn,">",cutoff))
  }
  
  
  xy <- cbind(pp$x,pp$y)
  n <- nrow(xy)
  Lam <- attr(pp, 'Lambda')
  summary(as.vector(rf.s <- log(Lam$v)))
  lgcp<-raster(as.im(list(x=Lam$xcol, y=Lam$yrow, z=t(exp(rf.s)))))
  plot(lgcp,axes=FALSE,box=FALSE,col=colo,legend.shrink=0.9,legend.width=4)
  #points(xy,pch=16,col=gray(0,0.01),cex=1)
  mtext(side=3,line=0,outer=FALSE,cex=0.6,text=paste("Simulated intensity",nrow(xy),"obs"))
  
  p<-st_as_sf(pp)
  st_crs(p)<-prj
  occ<-p[p$label!="window",]
  
  r<-raster(extent(win),ncol=100,nrow=100,crs=st_crs(p)$input)
  r<-extend(r,c(-20,1020,-20,1020))
  rx<-setValues(r,coordinates(r)[,1])
  ry<-setValues(r,coordinates(r)[,2])
  r<-stack(rx,ry)
  
  if(seed_effort){set.seed(seed)}
  model<-RMexp(var=params$evar[k],scale=params$erange[k])+RMnugget(var=0)+RMtrend(mean=params$emean[k])
  x.seq<-seq(extent(r)[1],extent(r)[2],length=300) 
  y.seq<-seq(extent(r)[3],extent(r)[4],length=300)
  sims<-RFsimulate(model,x=x.seq,y=y.seq)
  rsims<-raster(sims)#,xmn=extent(r)[1],xmx=extent(r)[2],ymn=extent(r)[3],ymx=extent(r)[4])
  modsims<-rsims
  modsims<-setValues(modsims,exp(params$eexp[k]*rescale(params$eb1[k]*coordinates(modsims)[,1]+params$eb2[k]*coordinates(modsims)[,2],1:0)))
  #par(mfrow=c(1,3))
  plot(modsims,axes=FALSE,box=FALSE,zlim=c(0,1),legend.shrink=0.9,legend.width=4,col=colorRampPalette(c("white","black"))(200))
  mtext(side=3,line=0,outer=FALSE,cex=0.6,text="Proportional effort gradient")
  eff<-exp(rsims)
  #plot(eff)
  eff<-exp(rsims)*modsims
  eff<-round(resample(eff,r,method="ngb"))
  eff2<-eff
  eff2[eff2==0]<-NA
  plot(eff2,axes=FALSE,box=FALSE,legend.shrink=0.9,legend.width=4,col=colorRampPalette(c("grey75","tomato","darkred","black"))(200),colNA="grey95")
  mtext(side=3,line=0,outer=FALSE,cex=0.6,text="Realized effort")
  #par(mfrow=c(1,1))
  eff<-stack(eff,setValues(eff,1:ncell(eff)))
  le<-raster::values(log(eff[[1]]))
  le<-ifelse(is.infinite(le),0,le)           
  r<-stack(r,eff[[1]],setValues(eff[[1]],le))
  names(r)<-c("xx","yy","effort","logeffort")
  
  e<-extract(eff,occ)
  occ$id<-e[,2]
  occ$effort<-e[,1]
  occ<-occ[order(occ$id),]
  
  if(seed_thinning){set.seed(seed)}
  a<-aggregate(label~id+effort,data=occ,length)
  a<-a[order(a$id),]
  l<-lapply(seq_len(nrow(a)),function(i){
    if(a$effort[i]==0){
      rep(FALSE,a$label[i])
    }else{
      pkeep<-a$effort[i]/max(occ$effort)
      k<-sample(0:1,a$label[i],prob=c(1-pkeep,pkeep),replace=TRUE)  
      as.logical(k)
    }
  })
  
  
  occ_thin<-occ[unlist(l),] # not 100% sure order is maintained
  #occ_thin<-occ # do not thin
  occ_thin<-st_crop(occ_thin,setNames(win[c(1,3,2,4)],c("xmin","ymin","xmax","ymax")))
  nrow(occ_thin)
  occ_thinsp<-as(occ_thin,"Spatial")
  plot(eff[[1]],col="grey95",box=FALSE,axes=FALSE,legend=TRUE,legend.shrink=0.9,legend.width=4)
  plot(occ_thinsp,pch=16,col=gray(0,0.15),add=TRUE)
  mtext(side=3,line=0,outer=FALSE,cex=0.6,text=paste("Thinned occurences",nrow(occ_thin),"obs"))
  
  region<-as(extent(win),"SpatialPolygons")
  proj4string(region)<-prj
  
  pedge<-0.01
  edge<-min(c(diff(bbox(region)[1,])*pedge,diff(bbox(region)[2,])*pedge))
  
  
  Mesh<-inla.mesh.2d(loc.domain = as(extent(region),"SpatialPoints"),max.edge = c(edge,edge*1.5),min.angle = 21,cutoff = edge/2,offset = c(edge,edge*2),crs = crs(region))
  #explana<-explanaMesh(sPoly=region,mesh=Mesh,X=r)
  explana<-explanaMesh(sPoly=st_as_sf(region),meshSpace=Mesh,meshTime=NULL,X=rast(r))
  weight<-ppWeight(sPoly=st_as_sf(region), mesh=Mesh)
  
  
  bpriors<-list(prec=list(default=1/(10)^2,Intercept=1/(100)^2,effort=1/(10)^2,logeffort=1/(0.00000001)^2),mean=list(default=0,Intercept=0,effort=0,logeffort=1))
  m<-ppSpace(y ~ xx+yy, sPoints = occ_thin,
             explanaMesh = explana,
             ppWeight = weight,
             prior.range = c(50,0.01),
             prior.sigma = c(1,0.01),
             num.threads = 7,
             many = TRUE,
             control.inla = list(int.strategy = "eb"),
             fix = NULL,
             sboffset = "effort",
             control.fixed = bpriors,
             control.compute=list(config=TRUE),
             orthoCons = FALSE,
             verbose=FALSE
  )
  
  
  ml[[k]]<-m
  xlim<-bbox(region)[1,]
  ylim<-bbox(region)[2,]
  #mapMean<-mapSpace(ml[[k]],dims=dim(r)[1:2]*1,type="mean")
  #mapMean <- mask(mapMean,region) #spacePoly
  #mapSd<-mapSpace(ml[[k]],dims=dim(r)[1:2]*1,type="sd")
  #mapSd <- mask(mapSd,region) #spacePoly
  maps<-mapSpace(ml[[k]],dims=dim(r)[1:2]*1,sample=FALSE)
  maps <- mask(maps,region) #spacePoly
  mapMean <- maps[["mean"]]
  mapSd <- maps[["sd"]]
  plot(mapMean, col = colo, axes = FALSE, box = FALSE,xlim=xlim,ylim=ylim,legend.shrink=0.9,legend.width=4)
  mtext(side=3,line=0,outer=FALSE,cex=0.6,text="Predicted intensity")
  plot(mapSd, col = colorRampPalette(c("grey90","tomato","darkred","black"))(200), axes = FALSE, box = FALSE,xlim=xlim,ylim=ylim,legend.shrink=0.9,legend.width=4)
  mtext(side=3,line=0,outer=FALSE,cex=0.6,text="Sd intensity")
  #plot(inla.mesh2sp(Mesh)$triangles,border=adjustcolor("white",0.1),add=TRUE)
  plot(attributes(weight)$dmesh,border=adjustcolor("white",0.1),add=TRUE)
  #grid(10,6,col=gray(0,0.25),lty=3)
  summary(m)$coef
  print(seed)
  
}

lapply(ml,function(i){ summary(i)$coef})
dev.off()
file.show("C:/Users/God/Downloads/simulated_effort.png")
