
library(mapSpecies)
library(sf)
library(spatstat)
library(fields)
library(scales)
library(raster)
library(RandomFields)

###############################################
### rLGCP
###############################################

cat("\014")
seed<-sample(1:1000,1)
set.seed(seed)
print(seed)

colo<-colorRampPalette(c("grey85","steelblue4","steelblue2","gold2","tomato2","red4"))(200)
prj<-"+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"


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
  rangex=c(400), # range of spatial field
  #kappa=c(sqrt(8*nu)/rangex), # matern smoothness parameter (2) and range of spatial field 
  b0=c(-5,-7,-5,-7), # intercept (determines the nb of points generated
  b1=c(0.002), # beta coefficient on X gradient
  b2=c(0.002), # beta coefficient on Y gradient
  sigma=sqrt(c(0.00000001,0.00000001,0.5,0.5)), # sd of spatial field
  
  ### Simulated effort
  evar=c(2), # variance of effort spatial field
  erange=c(25), # range of effort spatial field
  emean=c(-1), # mean of effort spatial field
  eb1=c(-1), # strength of gradient in x
  eb2=c(-1), # strength of gradient in y
  eexp=c(-1) # exponential decrease in effort gradient
)

params<-lapply(params,function(i){rep(i,length.out=max(sapply(params,length)))})

#layout(matrix(seq_len(5*length(params[[1]])),nrow=5,byrow=FALSE))

par(mfcol=c(length(params[[1]]),5),oma=c(0,4,0,3),mar=c(1,1,1,1))

for(k in seq_along(params[[1]])){
  
  if(seed_intensity){set.seed(seed)}
  pp<-rLGCP("matern",mu=function(x,y){params$b0[k]+params$b1[k]*x+params$b2[k]*y},var=params$sigma[k]^2,scale=1/c(sqrt(8  *params$nu[k])/params$rangex[k]),nu=params$nu[k],win=o,dimyx=c(300,300))
  cutoff<-300000
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
  rsims<-raster(sims,xmn=extent(r)[1],xmx=extent(r)[2],ymn=extent(r)[3],ymx=extent(r)[4])
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
  plot(eff2,axes=FALSE,box=FALSE,legend.shrink=0.9,legend.width=4,col=colorRampPalette(c("grey90","tomato","darkred","black"))(200))
  mtext(side=3,line=0,outer=FALSE,cex=0.6,text="Realized effort")
  #par(mfrow=c(1,1))
  eff<-stack(eff,setValues(eff,1:ncell(eff)))
  le<-values(log(eff[[1]]))
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
  occ_thin<-st_crop(occ_thin,setNames(win[c(1,3,2,4)],c("xmin","ymin","xmax","ymax")))
  nrow(occ_thin)
  occ_thinsp<-as(occ_thin,"Spatial")
  plot(eff[[1]],col="white",box=FALSE,axes=FALSE,legend=TRUE,legend.shrink=0.9,legend.width=4)
  plot(occ_thinsp,pch=16,col=gray(0,0.15),add=TRUE)
  mtext(side=3,line=0,outer=FALSE,cex=0.6,text=paste("Thinned occurences",nrow(occ_thin),"obs"))
  
  region<-as(extent(win),"SpatialPolygons")
  proj4string(region)<-prj
  
  pedge<-0.05
  edge<-min(c(diff(bbox(region)[0.2,])*pedge,diff(bbox(region)[2,])*pedge))
  
  
  Mesh<-inla.mesh.2d(loc.domain = as(extent(region),"SpatialPoints"),max.edge = c(edge,edge*1.5),min.angle = 21,cutoff = edge/2,offset = c(edge,edge*2),crs = crs(region))
  explana<-explanaMesh(sPoly=region,mesh=Mesh,X=r)
  weight<-ppWeight(sPoly=region, mesh=Mesh)
  
  
  bpriors<-list(prec=list(default=1/(10)^2,Intercept=1/(100)^2,effort=1/(10)^2,logeffort=1/(0.00000001)^2),mean=list(default=0,Intercept=0,effort=0,logeffort=1))
  m<-ppSpace(y ~ xx+yy, sPoints = occ_thinsp,
             explanaMesh = explana,
             ppWeight = weight,
             prior.range = c(100,0.05),
             prior.sigma = c(1,0.05),
             num.threads = 7,
             many = TRUE,
             control.inla = list(int.strategy = "eb"),
             fix = NULL,
             sboffset = "effort",
             control.fixed = bpriors,
             orthoCons = FALSE
  )
  
  
  m<-list(m=m)
  xlim<-bbox(region)[1,]
  ylim<-bbox(region)[2,]
  mapMean<-mapSpace(m[[1]],dims=dim(r)[1:2],type="mean")
  mapMean <- mask(mapMean,region) #spacePoly
  plot(mapMean, col = colo, axes = FALSE, box = FALSE,xlim=xlim,ylim=ylim,legend.shrink=0.9,legend.width=4)
  mtext(side=3,line=0,outer=FALSE,cex=0.6,text="Predicted intensity")
  #grid(10,6,col=gray(0,0.25),lty=3)
  lapply(m,function(i){summary(i)$coef})
  print(seed)
  
}

