

library(mapSpecies)
library(sf)
library(data.table)
library(scales)
library(raster)

########################
########################

r<-stack(predictors[[names(predictors)[names(predictors)%in%c("tmean","tmean2","broadleafs","broadleafs2","sbias")]]])


bpriors<-list(prec=list(default=1/(0.3)^2,latitude=1/(10)^2,latitude2=1/(10)^2,Intercept=1/(100)^2,sbias=1/(100)^2),mean=list(default=0,Intercept=0,latitude=0,latitude2=0,sbias=0))

pedge<-0.025
edge<-min(c(diff(bbox(region)[0.2,])*pedge,diff(bbox(region)[2,])*pedge))

Mesh <- inla.mesh.2d(loc.domain = coordinates(occsp),
                     max.edge = c(edge,edge*1.5),
                     min.angle = 21,
                     cutoff = edge/2,
                     offset = c(edge,edge*2),
                     crs = crs(region))

explana<-explanaMesh(sPoly=region,mesh=Mesh,X=r)

weight <- ppWeight(sPoly = region, mesh = Mesh)


mpp <- ppSpace(y ~ tmean+tmean2, sPoints = occsp,
               explanaMesh = explana,
               ppWeight = weight,
               #prior.range = c(200,0.01),
               #prior.sigma = c(0.1,0.01),
               num.threads=7,
               many=TRUE,
               control.inla = list(int.strategy = "eb"),
               fix = NULL,
               sboffset = "sbias",
               control.fixed = bpriors
)


m<-list(mpp=mpp)
par(mfrow=n2mfrow(length(m)),mar=c(3,3,3,7))
for(i in 1:length(m)){
  mapMean<-mapSpace(m[[i]],dims=dim(r)[1:2],type="mean",sPoly=as(q,"Spatial"))
  mapMean <- mask(mapMean,q) #spacePoly
  #mapMean<-crop(mapMean,extent(-350,1000,-300,800))
  names(mapMean)<-paste(sp,names(m)[i])
  init[[ninit]]<-mapMean
  ninit<-ninit+1
  plot(mapMean, col = colo, axes = TRUE, box = FALSE, main = paste(names(m)[i],sp,sep=" - "))
  #points(pressp,pch=16,col=alpha("forestgreen",0.99),cex=0.2)
  points(occsp,pch=16,col=alpha("black",0.05),cex=0.6)
  plot(st_geometry(q),add=TRUE,border=gray(0,0.3))
  ### contour
  ID<-inla.stack.index(attributes(m[[i]])$Stack,tag="est")$data
  #contour(mapMean,levels=quantile(m[[i]]$summary.fitted.values[["mean"]][ID],0.05),add=TRUE,col="red")
}
par(mfrow=c(1,1))


