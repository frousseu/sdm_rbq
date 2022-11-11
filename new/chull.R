
#a1<-st_read("/data/predictors_sdm/expert_maps/IUCN/REPTILES.shp")
#a2<-st_read("/data/predictors_sdm/expert_maps/IUCN/AMPHIBIANS.shp")
#a<-rbind(a1,a2)
#a$species<-a$binomial
#em<-st_transform(a,st_crs(s))

# https://www.usgs.gov/programs/gap-analysis-project/science/species-data-download

splist<-d$species
#splist<-splist[splist%in%em$species]
sp<-sample(names(rev(sort(table(splist))))[1:length(unique(splist))],1)
#sp<-"Lithobates palustris"
x<-s[s$species==sp,]
plot(st_geometry(x),pch=16,col=adjustcolor("darkgreen",0.2))
plot(st_geometry(na),add=TRUE,border=gray(0,0.15))

n<-1:nrow(x)
n<-n[unique(round(seq(min(n),max(n),length.out=100),0))]
l<-lapply(n,function(i){
  print(i)
  xx<-x[sample(1:nrow(x),i,replace=FALSE),]
  st_union(xx)
})

#obs<-rev(sort(table(x$recordedBy)))[1:10]

m<-do.call("c",l)
m<-st_convex_hull(m)
par(mfrow=c(2,2),mar=c(2,2,0,0))
plot(st_buffer(st_geometry(x),dist=250),border="white")
w<-which(em$species==sp)
if(any(w)){
  plot(st_geometry(em[w,]),add=TRUE,col=gray(0,0.15),border=NA)
}
plot(st_geometry(x),pch=16,col=adjustcolor("darkgreen",0.2),add=TRUE)
plot(st_geometry(na),border=gray(0,0.15),add=TRUE)
plot(m,border=gray(0,0.05),add=TRUE)
mtext(side=3,line=-1.5,adj=0.0,text=sp,font=2,xpd=TRUE)

a <- as.numeric(max(st_area(m)))
b <- 0.01
model <- nls(Y ~ ((a/b)*X)/(1+(X/b)),
 data = data.frame(Y=as.numeric(st_area(m)),X=n-1),
 start = list(a = a, b = b))
 plot(n,st_area(m),ylim=range(c(0,coef(model)["a"],st_area(m)))*1.05)
lines(0:max(n),predict(model,data.frame(X=0:max(n))),col="red",lwd=2)
abline(h=coef(model)["a"],lty=3,col="red")
mtext(side=1,adj=0.9,text=round(a/coef(model)["a"],2),font=2,col="red",line=-3,cex=2)

hist(st_area(m),breaks=50)


rangeHull<-function(x,species="species",breaks=50,plot=TRUE){
  n<-1:nrow(x)
  n<-n[unique(round(seq(min(n),max(n),length.out=breaks),0))]
  l<-lapply(n,function(i){
    print(i)
    xx<-x[sample(1:nrow(x),i,replace=FALSE),]
    st_union(xx)
  })

#obs<-rev(sort(table(x$recordedBy)))[1:10]

  m<-do.call("c",l)
  m<-st_convex_hull(m)
  
  par(mfrow=c(2,2),mar=c(2,2,0,0))
  plot(st_buffer(st_geometry(x),dist=250),border="white")
  #w<-which(em$species==i)
  #if(any(w)){
  #  plot(st_geometry(em[w,]),add=TRUE,col=gray(0,0.15),border=NA)
  #}
  plot(st_geometry(x),pch=16,col=adjustcolor("darkgreen",0.2),add=TRUE)
  plot(st_geometry(na),border=gray(0,0.15),add=TRUE)
  plot(m,border=gray(0,0.05),add=TRUE)
  mtext(side=3,line=-1.5,adj=0.0,text=species,font=2,xpd=TRUE)

  a <- as.numeric(max(st_area(m)))
  b <- 0.01
  model<-tryCatch(nls(Y ~ ((a/b)*X)/(1+(X/b)),
    data = data.frame(Y=as.numeric(st_area(m)),X=n-1),
    start = list(a = a, b = b)),error=function(j){TRUE})
  if(!isTRUE(model)){  
    plot(n,st_area(m),ylim=range(c(0,coef(model)["a"],st_area(m)))*1.05)
    lines(0:max(n),predict(model,data.frame(X=0:max(n))),col="red",lwd=2)
    abline(h=coef(model)["a"],lty=3,col="red")
    mtext(side=1,adj=0.9,text=round(a/coef(model)["a"],2),font=2,col="red",line=-3,cex=2)
    hist(st_area(m),breaks=50)
    list(hull=m,nobs=nrow(x),a=coef(model)["a"],reach=a/coef(model)["a"])
  }else{
    list(hull=m,nobs=nrow(x),a=NA,reach=NA)
  }
}

rangeHull(loccs[[1]],species=names(loccs)[1],breaks=50)


tab<-table(d$species)
tab<-tab[tab>50]
sp<-sample(names(tab),100)
r<-sapply(sp,function(i){
  print(i)
  #splist<-splist[splist%in%em$species]
  #sp<-sample(names(rev(sort(table(splist))))[1:length(unique(splist))],1)
#sp<-"Lithobates palustris"
  x<-s[s$species==i,]
  r<-rangeHull(x,breaks=20)
  r$reach
})
  
prop<-length(which(r>=0.95))/length(r) 
prop*length(tab)
prop*length(tab)/length(unique(d$species))



1+1