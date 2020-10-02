
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

seed<-sample(1:1000,1)
set.seed(seed)

colo<-colorRampPalette(c("grey92","steelblue4","steelblue2","gold2","tomato2","red4"))(200)
prj<-"+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

win<-c(0,1000,0,1000)
o<-owin(win[1:2],win[3:4])
nu=2
rangex=300
kappa=sqrt(8*nu)/rangex
b0=-6
b1=0.002
b2=b1
sigma=sqrt(0.00000001)
#exp(b0) * diff(range(o$x)) * diff(range(o$y))

#set.seed(1234) # 123 avec -17 aberrant
p<-rLGCP("matern",mu=function(x,y){b0+b1*x+b2*y},var=sigma^2,scale=1/kappa,nu=nu,win=o,dimyx=c(300,300))
cutoff<-100000
if(p$n>cutoff){
  pn<-p$n
  rm(p)
  stop(paste("Too many points generated",pn,">",cutoff))
}
plot(p,axes=FALSE);p

xy <- cbind(p$x,p$y)
n <- nrow(xy)
Lam <- attr(p, 'Lambda')
summary(as.vector(rf.s <- log(Lam$v)))
par(mfrow=c(1,2),oma=c(0,4,0,4))
image.plot(list(x=Lam$xcol, y=Lam$yrow, z=t(exp(rf.s))), main='log-Lambda',asp=1,col=colo)
points(xy,pch=16,col=gray(0,0.15),cex=0.6)

p<-st_as_sf(p)
st_crs(p)<-prj#32617
#p<-st_transform(p,prj)
occ<-p[p$label!="window",]

r<-raster(extent(win),ncol=100,nrow=100,crs=st_crs(p)$input)
r<-extend(r,c(-200,1200,-200,1200))
rx<-setValues(r,coordinates(r)[,1])
ry<-setValues(r,coordinates(r)[,2])
r<-stack(rx,ry)

model<-RMexp(var=1.5,scale=25)+RMnugget(var=0)+RMtrend(mean=0)
#x.seq<-seq(o$xrange[1],o$xrange[2],length=300) 
#y.seq<-seq(o$yrange[1],o$yrange[2],length=300)

x.seq<-seq(extent(r)[1],extent(r)[2],length=300) 
y.seq<-seq(extent(r)[3],extent(r)[4],length=300)
sims<-RFsimulate(model,x=x.seq,y=y.seq)
eff<-exp(raster(sims))
#hist(values(eff))
eff<-round(resample(eff,r,method="ngb"))
plot(eff)
eff<-stack(eff,setValues(eff,1:ncell(eff)))
le<-values(log(eff[[1]]))
le<-ifelse(is.infinite(le),-1,le)           
r<-stack(r,eff[[1]],setValues(eff[[1]],le))
names(r)<-c("xx","yy","effort","logeffort")

e<-extract(eff,occ)
occ$id<-e[,2]
occ$effort<-e[,1]
occ<-occ[order(occ$id),]

l<-split(occ,occ$id)
l<-lapply(l,function(i){
  pkeep<-i$effort[1]/max(occ$effort)
  k<-sample(0:1,nrow(i),prob=c(1-pkeep,pkeep),replace=TRUE)  
  #i[as.logical(k),]
  as.logical(k)
})
#occ_thin<-do.call("rbind",l)
occ_thin<-occ[unlist(l),] # not 100% sure order is maintained
#occ_thin<-st_crop(occ_thin,c(xmin=0,ymin=0,xmax=1000,ymax=1000))
occ_thin<-st_crop(occ_thin,setNames(win[c(1,3,2,4)],c("xmin","ymin","xmax","ymax")))
nrow(occ_thin)
occ_thinsp<-as(occ_thin,"Spatial")

region<-as(extent(win),"SpatialPolygons")
proj4string(region)<-prj

pedge<-0.04
edge<-min(c(diff(bbox(region)[0.2,])*pedge,diff(bbox(region)[2,])*pedge))
edge

bpriors<-list(prec=list(default=1/(10)^2,latitude=1/(10)^2,latitude2=1/(10)^2,Intercept=1/(100)^2,sbias=1/(100)^2,effort=1/(100)^2,logeffort=1/(100)^2),mean=list(default=0,Intercept=0,latitude=0,latitude2=0,sbias=0,effort=0,logeffort=0))

Mesh<-inla.mesh.2d(loc.domain = as(extent(region),"SpatialPoints"),max.edge = c(edge,edge*1.5),min.angle = 21,cutoff = edge/2,offset = c(edge,edge*2),crs = crs(region))
explana<-explanaMesh(sPoly=region,mesh=Mesh,X=r)
weight<-ppWeight(sPoly=region, mesh=Mesh)

m1<-ppSpace(y ~ xx+yy+effort, sPoints = occ_thinsp,
            explanaMesh = explana,
            ppWeight = weight,
            prior.range = c(20000,NA),
            prior.sigma = c(0.00001,NA),
            num.threads = 7,
            many = TRUE,
            control.inla = list(int.strategy = "eb"),
            bias = "effort",
            control.fixed = bpriors#,#,#,
            #orthoCons = TRUE
)


m2<-ppSpace(y ~ xx+yy+offset(logeffort), sPoints = occ_thinsp,
             explanaMesh = explana,
             ppWeight = weight,
             prior.range = c(200,0.5),
             prior.sigma = c(1,0.5),
             num.threads = 7,
             many = TRUE,
             control.inla = list(int.strategy = "eb"),
             bias = NULL,
             control.fixed = bpriors,
             orthoCons = TRUE
)


m<-list(m1=m1,m2=m2)
par(mfrow=c(2,2),oma=c(0,4,0,4))
xlim<-bbox(region)[1,]
ylim<-bbox(region)[2,]
image.plot(list(x=Lam$xcol, y=Lam$yrow, z=t(exp(rf.s))), main=paste0('Simulated intensity - ',nrow(occ)," obs",sep=" "),asp=1,col=colo)
#points(xy,pch=1,col=gray(0,0.05),cex=0.6)
grid(10,6,col=gray(0,0.25),lty=3)
plot(eff[[1]],xlim=xlim,ylim=ylim,main="Effort (en entier, beaucoup de 0)",col=colo)
grid(10,6,col=gray(0,0.25),lty=3)
mapMean<-mapSpace(m[[1]],dims=dim(r)[1:2],type="mean")
mapMean <- mask(mapMean,region) #spacePoly
plot(mapMean, col = colo, axes = TRUE, box = FALSE, main = paste0(names(m)[1]," - ",nrow(occ_thin)," obs",sep=" "),xlim=xlim,ylim=ylim,legend.shrink=0.9,legend.width=1.5)
points(occ_thinsp,pch=1,col=gray(0,0.25),cex=0.6)
grid(10,6,col=gray(0,0.25),lty=3)
mapMean<-mapSpace(m[[2]],dims=dim(r)[1:2],type="space")
mapMean <- mask(mapMean,region) #spacePoly
plot(mapMean, col = colo, axes = TRUE, box = FALSE, main = paste0(names(m)[2]," - ",nrow(occ_thin)," obs",sep=" "),xlim=xlim,ylim=ylim,legend.shrink=0.9,legend.width=1.5)
points(occ_thinsp,pch=1,col=gray(0,0.25),cex=0.6)
grid(10,6,col=gray(0,0.25),lty=3)
lapply(m,function(i){summary(i)$coef})
print(seed)



#ids<-inla.stack.index(attributes(mppSpace)$Stack, tag="pred")$data
#vs<-mppSpace$summary.fitted.values[["mean"]][ids]
#rv<-setValues(r[[1]],vs)
#plot(rv)

#rr<-rasterize(occ,r[[1]],fun="count",background=0)
rr<-rasterize(occ_thin[,1],r[[1]],fun="count",background=0)


temp<-data.frame(eff=values(eff[[1]]),logeff=log(values(eff[[1]])),cou=values(rr)[,1],xx=coordinates(rr)[,1],yy=coordinates(rr)[,2])
temp2<-temp[temp$eff>0,]
fit<-glm(cou~xx+yy+logeff,data=temp2,family="poisson")
summary(fit)
newdat<-temp
newdat$eff<-max(temp$eff)
p<-predict(fit,newdat,type="response")
rrp<-setValues(rr[[1]],p)
plot(rrp,col=colo)
newdat<-data.frame(xx=500,yy=500,logeff=seq(0,max(temp$logeff),by=0.01))
p<-predict(fit,newdat,type="response")
plot(jitter(temp$logeff),jitter(temp$cou))
lines(newdat$logeff,p)



temp<-data.frame(eff=values(eff[[1]]),logeff=log(values(eff[[1]])),cou=values(rr)[,1],xx=coordinates(rr)[,1],yy=coordinates(rr)[,2])
temp2<-temp[temp$eff>0,]
fit<-glm(cou~xx+yy+log(eff),data=temp2,family="poisson")
summary(fit)
#visreg(fit,"logeff")
#visreg(fit,"logeff",scale="response")
plot(jitter(temp2$eff),jitter(temp2$cou),xlim=c(0,max(temp2$eff)))
newdat<-data.frame(xx=500,yy=500,eff=seq(0,max(temp2$eff),by=0.01))
p<-predict(fit,newdat,type="response")
lines(newdat$eff,p)

sims<-simulateResiduals(fit)
plot(sims)



rrp2<-setValues(rrp,rescale(values(rrp),0:1)
plot(rrp2-rrp1)
grid()
                
plot(values(eff[[1]]),values(rr)[,1])
temp<-data.frame(eff=values(eff[[1]]),cou=values(rr)[,1])
fit<-lm(cou~eff,data=temp)
plot(values(eff[[1]]),values(rr)[,1])
abline(fit)
                
                
                