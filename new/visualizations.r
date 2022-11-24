
library(tidyr)
library(foreach)
library(doParallel)


options(vsc.dev.args = list(width = 1200, height = 800))

#############################################################
### Plot effects ############################################

#plot(exp(log(sdm[["mean"]])-sdm[["spacemean"]]))
 plot(exp(sdm$linkmean-sdm$spacemean))

#listo<-lapply(ls(), function(x) object.size(get(x)))
#o<-rev(order(unlist(listo)))
#setNames(sapply(listo,format,units="auto")[o],ls()[o])

#############################################################
### Marginal effects  #######################################

checkpoint("Marginal effects")

m<-mpp
class(m)<-"inla"
nsims<-500
samples<-inla.posterior.sample(nsims,m)

params<-dimnames(m$model.matrix)[[2]]
nparams<-sapply(params,function(i){
  match(paste0(i,":1"),row.names(samples[[1]]$latent)) 
}) 
tab<-table(sapply(strsplit(row.names(samples[[1]]$latent),":"),"[",1))
tab
nweights<-grep("i:",row.names(samples[[1]]$latent))

nd<-newdata2(as.data.frame(as.matrix(m$model.matrix)[,-1,drop=FALSE])[1:(tab[1]/2),,drop=FALSE],n=100,n2=5,fun=mean)
v<-names(nd)
nd<-unlist(nd,recursive=FALSE)

marginals<-lapply(seq_along(nd),function(k){
  p<-lapply(1:nsims,function(i){
    betas<-samples[[i]]$latent[nparams]
    names(betas)<-names(nparams)
    fixed<-cbind(Intercept=1,as.matrix(nd[[k]]))[,names(betas)] %*% betas # make sure betas and vars are in the same order
    # this if we want a spatial part
    #wk<-samples[[i]]$latent[nweights]
    #if(is.factor(xs@data[,v[k]])){ # factors never in model (et)
    #spatial<-as.matrix(inla.spde.make.A(mesh=mesh,loc=matrix(c(0.3,0.5),ncol=2)[rep(1,nlevels(size[,v[k]])),,drop=FALSE])) %*% wk
    #}else{
    #spatial<-as.matrix(AA) %*% wk # stack was Apn in fire
    #}
    #p<-fixed # ignores spatial part
    list(fixed,fixed)
  })
  p1<-do.call("cbind",lapply(p,"[[",1))
  p2<-do.call("cbind",lapply(p,"[[",2))
  p1<-t(apply(p1,1,function(i){c(quantile(i,0.025),mean(i),quantile(i,0.975))}))
  p2<-t(apply(p2,1,function(i){c(quantile(i,0.025),mean(i),quantile(i,0.975))}))
  p1<-exp(p1)
  p2<-exp(p2)
  p1<-setNames(as.data.frame(p1),c("low","mean","high"))
  add<-strsplit(names(nd)[k],"\\.")[[1]]
  data.frame(vars=add[1],block=add[2],p1)
})
names(marginals)<-names(nd)


preds<-cbind(do.call("rbind",marginals),do.call("rbind",nd))
ma<-match(names(preds),row.names(means))
ma<-ma[!is.na(ma)]
preds[row.names(means)[ma]]<-lapply(row.names(means)[ma],function(i){
  backTransform(preds[[i]],i)
})


preds<-split(preds,preds$vars)
preds<-lapply(preds,function(i){
  if(all(is.na(i$block))){
    list(i)
  }else{
    split(i,i$block)
  }
})


png(file.path(getwd(),"graphics","marginal_effects.png"),width=10,height=8,units="in",res=300)
par(mfrow=n2mfrow(length(preds),asp=2),mar=c(3,2,1,1),oma=c(1,1,1,1))
lapply(names(preds),function(i){
  xvar<-strsplit(i,"_")[[1]][1]
  xvar2<-strsplit(i,"_")[[1]][2]
  xlim<-range(preds[[i]][[1]][,xvar])
  ywide<-unlist(preds)
  ylim<-range(as.numeric(ywide[grep("\\.mean",names(ywide))]))
  #ylim<-range(unlist(do.call("rbind",preds[[i]])[,c("low","mean","high"),drop=FALSE]))
  plot(0.1,0.1,xlab="",ylab="",xlim=xlim,ylim=ylim,type="n",mgp=c(2,0.5,0),tcl=-0.1,log="")
  ### density of values where species is observed
  vals<-backTransform(dmeshPred[,i],i)
  h<-hist(vals[dmesh$spobs>0],breaks=seq(min(vals),max(vals),length.out=40),plot=FALSE)
  h$density<-(h$density/max(h$density))
  h$density<-h$density/(max(h$density)/ylim[2])
  lapply(seq_along(h$density),function(j){
    rect(xleft=h$breaks[j],ybottom=0,xright=h$breaks[j+1],ytop=h$density[j],border=NA,col=adjustcolor("seagreen",0.25))
  })

  lapply(seq_along(preds[[i]]),function(jj){
    j<-preds[[i]][[jj]]
    x<-j[,xvar]
    polygon(c(x,rev(x),x[1]),c(j[,"low"],rev(j[,"high"]),j[,"low"][1]),col=adjustcolor("black",0.25),border=NA)
    lines(x,j[,"mean"],lwd=1.5,col=adjustcolor("black",0.25))
    mtext(side=1,line=2,text=xvar,font=2,cex=1.25)
  })
  if(length(preds[[i]])>1){
    temp<-do.call("rbind",preds[[i]])
    var2<-strsplit(temp$vars[1],"_")[[1]][2]
    legend("top",title=var2,legend=round(as.numeric(unique(temp[,var2])),2),lwd=1,col=seq_along(preds[[i]]),bty="n",cex=1.25)
  }





  #p1<-marginals[[i]]
  #ylim<-c(0,max(unlist(marginals)))
  #ylim<-c(0,max(as.vector(p1[,2])))
  #plot(vals,p1[,2],type="l",ylim=ylim,xlab=i,ylab="Intensity",font=1,lty=1,mgp=c(1.5,0.45,0),tcl=-0.3)
  #lines(vals,p1[,2],lwd=1.5)
  #polygon(c(vals,rev(vals),vals[1]),c(p1[,1],rev(p1[,3]),p1[,1][1]),col=alpha("black",0.1),border=NA)
})
dev.off()

#preds1<-mapSpace(ml[[1]],dims=dim(r),sPoly=region)[["mean"]]

##############
##############
##############
##############

### check vif for specific model
library(car)
vs<-unique(gsub("[[:digit:]]","",vars_pool))
vs<-vs[!vs%in%c("sbias","fixed")]
#vs<-c(vs,"water")
dat<-as.data.frame(dmeshPred)
dat$y<-rnorm(nrow(dat))
mo<-lm(formula(paste("y~",paste(vs,collapse="+"))),data=dat)
vif(mo)



##############
##############
##############




ebird<-rast("/data/predictors_sdm/expert_maps/eBird/abundance/Vireo_gilvus_ebird.tif")
ebird<-project(ebird,sdm)


dmesh$pred<-mpp$summary.linear.predictor[1:nrow(dmesh),1]
dmesh$sd<-mpp$summary.linear.predictor[1:nrow(dmesh),2]
dmesh$spacepred<-mpp$summary.random$i[1:nrow(dmesh),2]
plot(dmesh["spacepred"],pal=magma,nbreaks=100,lwd=0.05,border=adjustcolor("white",0.25),mar=c(0,0,0,0),key.pos=4,reset=FALSE)
plot(st_geometry(na),lwd=0.30,border=adjustcolor("white",0.85),mar=c(0,0,0,0),add=TRUE)



sdm<-mapSpace(modelSpace=m[[i]],dims=round(0.25*c(2,2)*dim(r)[1:2],0),sPoly=as(na,"Spatial"),sample=TRUE)#[[type]]

sdm2<-sdm[[grep("sample",names(sdm))]]
#sdm2<-aggregate(sdm2,5)
a<-as.array(exp(sdm2))
pixels<-sapply(unlist(asplit(asplit(a,1:2),1),recursive=FALSE,use.names=FALSE),function(i){
  #j<-i[i>=quantile(i,0.025) & i<=quantile(i,0.975)]
  sample(i,size=1)
})
pixelation<-setValues(sdm2[[1]],pixels)
pixelation<-resample(pixelation,sdm[[1]])
plot(mask(pixelation,vect(na)))


un<-list(
  sdm$mean,
  sdm$sd,
  sdm$linksd,
  sdm$spacesd,
  exp(sdm[["sample001"]]),
  pixelation

)
un<-rast(un)
names(un)<-c("abundance","sd on abundance","sd on the link scale","sd of the spatial component","a sample from the posterior abundance","a sample from the posterior of each pixel")
un<-mask(un,vect(na))

png(file.path(getwd(),"graphics","uncertainties.png"),width=10,height=6,units="in",res=400)
par(mfrow=n2mfrow(nlyr(un),asp=2),oma=c(0,0,0,0))
  lapply(seq_len(nlyr(un)),function(i){
  plot(un[[i]],col=rev(magma(200)),axes=FALSE,mar=c(0.5,0,1.5,4))
  plot(st_geometry(na),border=adjustcolor("black",0.25),lwd=0.51,add=TRUE)
  plot(st_geometry(nalakes),col="white",border=adjustcolor("black",0.15),lwd=0.51,add=TRUE)
  mtext(side=3,text=names(un)[i],font=2)
})
dev.off()

v<-values(sdm$spacemean)
v<-seq(min(v,na.rm=TRUE),max(v,na.rm=TRUE),length.out=200)
plot(sdm$spacemean,col=colo.scale(v,cols=c("dodgerblue4","dodgerblue","lightgoldenrod","firebrick3","firebrick4"),center=TRUE))


############################################
### show eBird
#list.files(tools::R_user_dir("ebirdst"),recursive=TRUE)
ebirdpath<-list.files("/data/predictors_sdm/expert_maps/eBird/abundance",full=TRUE,pattern="_ebird.tif")
ebird<-rast(ebirdpath[grep(gsub(" ","_",species),ebirdpath)])
#ebird<-project(ebird,crs(na))
#ebird<-crop(ebird,na)
#ebird<-sdm$mean
plot((ebird)^(1/1),col=colo[["mean"]])
plot(st_geometry(nalakes),border=gray(0,0.25),lwd=0.15,add=TRUE)
plot(st_geometry(na),border=gray(1,0.99),lwd=0.5,add=TRUE)
plot(st_geometry(st_union(na)),border=gray(0,0.15),lwd=0.5,add=TRUE)


##########################################################
##########################################################
##########################################################
### eBird comparison #####################################
##########################################################
##########################################################
##########################################################

lsdms<-list.files("/data/sdm_rbq/rasters",pattern="_birds.tif",full=TRUE)
lsdms<-lsdms[rev(order(file.info(lsdms)$mtime))]

df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
sp<-gsub(" ","_",df$species[which(df$reach>=0.00)])
lsdms<-lsdms[gsub("_birds.tif","",basename(lsdms))%in%sp]

lsdms<-lsdms[1:min(length(lsdms),5000)]
sdms<-lapply(lsdms,rast,lyrs="mean")
names(sdms)<-gsub("_birds.tif","",basename(lsdms))
dims<-sapply(sdms,function(i){paste(dim(i)[1:2],collapse="_")})
keep<-names(rev(sort(table(dims))))[1]
sdms<-rast(sdms[names(sdms)[dims%in%keep]])

matches<-match(gsub("_"," ",names(sdms)),d$species)
ebirdnames<-d$ebird[matches]
remove<-names(sdms)[is.na(ebirdnames)]
sdms<-sdms[[!names(sdms)%in%remove]] # remove what is not in ebird abundance 
ebirdnames<-d$ebird[match(gsub("_"," ",names(sdms)),d$species)]
ebirdnames<-gsub(" ","_",ebirdnames)

ebirdpath<-list.files("/data/predictors_sdm/expert_maps/eBird/abundance",full=TRUE,pattern="_ebird2.tif")
ebirdpath<-ebirdpath[unlist(sapply(ebirdnames,function(i){grep(i,ebirdpath)}),use.names=FALSE)]

ebird<-rast(lapply(ebirdpath,rast))
names(ebird)<-gsub("_ebird2.tif","",basename(ebirdpath))
#ebird<-project(ebird,sdms[[1]])
sdms<-sdms[[names(ebird)]]
#ebird<-project(ebird,crs(na))
#ebird<-crop(ebird,na)
#plot(ebird)

cors<-sapply(seq_len(nlyr(sdms)),function(i){
  print(i)
  #vsdms<-values(sdms[[i]])[,1]
  #vebird<-values(ebird[[i]])[,1]
  #cor(vsdms,vebird,use="complete.obs",method="pearson")
  #plot(vsdms,vebird)
  cor(values(sdms[[i]])[,1],values(ebird[[i]])[,1],use="complete.obs",method="pearson")
  #nicheOverlap(raster(sdms[[i]]),raster(ebird[[i]]),stat="I")
})
cors

### nicheOverlap
cl<-makeCluster(21)
registerDoParallel(cl)
no<-foreach(i=names(sdms),.packages=c("dismo","raster")) %dopar% {
  r1<-raster(file.path("/data/sdm_rbq/rasters",paste0(i,"_birds.tif")),layer=1)
  r2<-raster(file.path("/data/predictors_sdm/expert_maps/eBird/abundance",paste0(i,"_ebird2.tif")),layer=1)
  nicheOverlap(r1,r2,stat="I")
}
stopCluster(cl)
no<-unlist(no)


### adds latest correlation to df of results
res<-data.frame(
  species=gsub("_"," ",names(sdms)),
  correlation=cors,
  I=no
)
df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
m<-match(df$species,res$species)
df$correlation<-ifelse(!is.na(m),res$correlation[m],df$correlation)
df$I<-ifelse(!is.na(m),res$I[m],df$I)
write.csv(df,file="/data/sdm_rbq/graphics/mapSpeciesres.csv",row.names=FALSE,append=FALSE)

### plot it
png("/data/sdm_rbq/graphics/overlap.png",width=12,height=10,units="in",res=400)
ag<-df$I
ag<-ag[which(df$reach>=0.85)]
brks<-seq(0.0,1,by=0.1)
h<-hist(ag,xlim=c(0,1),breaks=brks,border="white",col=adjustcolor("seagreen",0.5),main="",xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,at=brks,mgp=c(2,1.25,0),pos=c(0,0),col="grey50",col.axis="black",col.ticks="grey50",lwd=5,cex.axis=1.5,font.axis=2)
axis(2,at=pretty(h$counts),mgp=c(2,1,0),pos=c(-0.0,0),col="grey50",col.axis="black",col.ticks="grey50",lwd=5,cex.axis=1.5,font.axis=2,las=2)
lines(c(rep(mean(ag,na.rm=TRUE),2)),c(0,par("usr")[4]),lwd=8,col="black")
mtext(side=1,line=2,text="Overlap",font=2,cex=3)
mtext(side=2,line=0.5,text="Number of species",font=2,cex=3)
dev.off()

names(sdms)[cors>-0.9 & cors<0.1]
i<-which(names(sdms)=="Anthus_rubescens")


### 
#sp<-"Pipilo_erythrophthalmus"
#r1<-raster(grep(sp,list.files("/data/sdm_rbq/rasters",pattern="\\.tif",full=TRUE),value=TRUE))
#r2<-raster(grep(sp,list.files("/data/predictors_sdm/expert_maps/eBird/abundance",pattern="_ebird2.tif",full=TRUE),value=TRUE))
#o<-nicheOverlap(r1,r1*1)
#o
#cor(values(r1),values(r2),use="complete.obs")

### used for using exactly the same area for all species
#common<-mean(sdms)+mean(ebird)
#sdms<-mask(sdms,common)
#sdms<-sdms/global(sdms,fun="max",na.rm=TRUE)[,1]
#ebird<-mask(ebird,common)
#ebird<-ebird/global(ebird,fun="max",na.rm=TRUE)[,1]

##################################################################
### Combinations
common<-raster(common)

splist<-rep(names(sdms)[1:10],2)
splist<-paste0(splist,rep(c("_S","_E"),each=length(splist)/2))
combs<-expand.grid(V1=splist,V2=splist,stringsAsFactors=FALSE)


### nicheOverlap
cl<-makeCluster(20)
registerDoParallel(cl)
no<-foreach(i=1:nrow(combs),.packages=c("dismo","raster")) %dopar% {
  if(any(grep("_S",combs[i,1]))){
    path1<-file.path("/data/sdm_rbq/rasters",gsub("_S","_birds.tif",combs[i,1]))
  }else{
    path1<-file.path("/data/predictors_sdm/expert_maps/eBird/abundance",gsub("_E","_ebird2.tif",combs[i,1]))
  }
  if(any(grep("_S",combs[i,2]))){
    path2<-file.path("/data/sdm_rbq/rasters",gsub("_S","_birds.tif",combs[i,2]))
  }else{
    path2<-file.path("/data/predictors_sdm/expert_maps/eBird/abundance",gsub("_E","_ebird2.tif",combs[i,2]))
  }
  r1<-raster(path1,layer=1)
  r2<-raster(path2,layer=1)
  nicheOverlap(mask(r1,common),mask(r2,common),stat="I")
}
stopCluster(cl)
no<-unlist(no)
combs$dist<-no
dat<-as.data.table(combs)
dat<-dat[order(V1,-dist),]


plan(multisession,workers=4)
cores<-nbrOfWorkers() # get nbr of workers from the chosen plan
chunks <- split(1:nrow(combs), rep(1:cores, each=ceiling(nrow(combs)/cores))[1:nrow(combs)])
options(future.globals.maxSize = 10000 * 1024 ^ 2)
combs$dist<-do.call("c",future_lapply(chunks,function(chunksi){
  sapply(chunksi,function(i){
    print(i)
    ZeroOneEuclideanDist(mat[[combs[i,1]]],mat[[combs[i,2]]])
  })
}))
plan(sequential)
#combs<-combs[combs$dist>0,]
dat<-as.data.table(combs)
dat$V1<-as.character(dat$V1)
dat$V2<-as.character(dat$V2)
dat<-dat[substr(dat$V1,nchar(dat$V1)-1,nchar(dat$V1))=="_s" & substr(dat$V2,nchar(dat$V2)-1,nchar(dat$V2))=="_e",]
dat<-dat[order(V1,dist),]








HellingerDist <- function (mat1,mat2)
 {
 p1 <- sum(mat1,na.rm=T)
 p2 <- sum(mat2,na.rm=T)
 return(sqrt(0.5*sum((sqrt(mat1/p1) - sqrt(mat2/p2))^2,na.rm=T)))
 }

NormEuclideanDist <- function(mat1,mat2)
 {
 p1 <- mat1/sum(mat1,na.rm = T)
 p2 <- mat2/sum(mat2,na.rm = T)
 return(sqrt(sum((p1 - p2)^2,na.rm=T)))
 }


ZeroOneEuclideanDist <- function(mat1,mat2)
 {
 p1 <- mat1/max(mat1,na.rm = T)
 p2 <- mat2/max(mat2,na.rm = T)
 return(sqrt(sum((p1 - p2)^2,na.rm=T)))
 }


mat<-lapply(c(sdms,ebird),as.matrix)
names(mat)<-c(paste0(names(sdms),"_s"),paste0(names(ebird),"_e"))

#combs<-as.data.frame(rbind(t(combn(names(mat),2)),matrix(rep(names(mat),2),ncol=2)))
combs<-expand.grid(V1=names(mat),V2=names(mat),stringsAsFactors=FALSE)

plan(multisession,workers=4)
cores<-nbrOfWorkers() # get nbr of workers from the chosen plan
chunks <- split(1:nrow(combs), rep(1:cores, each=ceiling(nrow(combs)/cores))[1:nrow(combs)])
options(future.globals.maxSize = 10000 * 1024 ^ 2)
combs$dist<-do.call("c",future_lapply(chunks,function(chunksi){
  sapply(chunksi,function(i){
    print(i)
    ZeroOneEuclideanDist(mat[[combs[i,1]]],mat[[combs[i,2]]])
  })
}))
plan(sequential)
#combs<-combs[combs$dist>0,]
dat<-as.data.table(combs)
dat$V1<-as.character(dat$V1)
dat$V2<-as.character(dat$V2)
dat<-dat[substr(dat$V1,nchar(dat$V1)-1,nchar(dat$V1))=="_s" & substr(dat$V2,nchar(dat$V2)-1,nchar(dat$V2))=="_e",]
dat<-dat[order(V1,dist),]


l<-split(dat,dat$V1)
names(l)<-gsub("_s","",names(l))
spe<-"Troglodytes_aedon"
g<-grep(spe,l[[spe]]$V2)
plot(c(sdms[[spe]],identity(ebird[[gsub("_e","",l[[spe]]$V2[1:g])]])))





spe1<-"Accipiter_cooperii"
spe2<-spe1
spe2<-"Hirundo_rustica"
mat1<-mat[[paste0(spe1,"_s")]]
mat2<-mat[[paste0(spe2,"_e")]]
mat1<-mat1/max(mat1,na.rm=TRUE)
mat2<-mat2/max(mat2,na.rm=TRUE)
HellingerDist(mat1,mat2)



distm <- dat %>% 
  spread(V2, dist, fill="") %>% 
  tibble::column_to_rownames(var="V1") %>% 
  as.matrix

m<-matrix(as.numeric(distm),ncol=ncol(distm))
colnames(m)<-colnames(distm)
rownames(m)<-rownames(distm)
#m<-m[,-1]
#m<-m[-nrow(m),]
m<-m[,order(colnames(m))]
m<-m[order(rownames(m)),]
m<-as.dist(m)


ord<-cmdscale(m,eig=TRUE,k=5)
plot(ord$points[,1:2],asp=1)
text(ord$points[,1:2],label=rownames(ord$points))
#barplot(ord$eig)

closest<-ord$points[substring(rownames(ord$points),nchar(rownames(ord$points))-1)=="_e",]
target<-ord$points[substring(rownames(ord$points),nchar(rownames(ord$points))-1)=="_s",]

nn<-knnx.index(closest,target,k=1)#[,1,drop=FALSE]

nn<-cbind(target=rownames(target),matrix(rownames(closest)[as.vector(nn)],ncol=ncol(nn))) |> as.data.frame()
nn

plot(ebird[[gsub("_e","",names(rev(sort(table(nn[,2])))))]])

apply(vals,2,function(i){table(is.na(i))})

plot(mean(ebird))
plot(st_geometry(na),add=TRUE)

mat1<-matrix(0:24,ncol=5)
mat2<-matrix(0:24,ncol=5)
mat2[1:15]<-NA
mat1[1:15]<-NA
mat1<-mat1/max(mat1,na.rm=TRUE)
mat2<-mat2/max(mat2,na.rm=TRUE)
HellingerDist(mat1,mat2)

sdms<-mask(sdms,mean(ebird))

maxs<-global(sdms,fun="max",na.rm=TRUE)


plot(sdms/maxs[,1])


plot(st_geometry(na))
plot(st_geometry(loccs[["Cathartes aura"]]),add=TRUE)


##### find names
library(rgbif)
library(taxize)
temp <- as.data.table(gnr_resolve(c("Circus cyaneus")))
head(temp)

x<-as.data.table(name_backbone(name="Gallinago wilsoni", rank='species'))
head(x)


spnames<-unique(d$species[!d$species%in%ed$species])

l<-lapply(spnames,function(i){
  #x<-name_backbone(name=i, rank='species')
  print(i)
  x <- gnr_resolve(i)
  #if(nrow(x)==0L){
  #  NULL
  #}else{
    res<-x[x$data_source_title==c("The eBird/Clements Checklist of Birds of the World"),]
    #res2<-x[x$data_source_title==c("Birds of the World: Recommended English Names"),]
    #res<-rbind(res1,res2)
    if(nrow(res)>0){
      if(i!=res$matched_name[1]){
        sp<-res$matched_name[1]
      }else{
        sp<-NA
      }
    }else{
      key<-name_backbone(i)$usageKey
      x<-name_usage(key=key,data="synonyms")$data$canonicalName
      k<-x%in%ed$species
      if(any(k)){
        sp<-x[k]
      }else{
        sp<-NA
      }
    }
    c(species=i,species_match=sp)
})
x<-as.data.table(do.call("rbind",l))

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
 

  #print(length(x))
  #x[x%in%ed$species]

})






x$species[is.na(x$species_match)]

ma<-match(x$species[is.na(x$species_match)],d$species)
d[ma,c("species","totobs")][order(-totobs),][1:50,]



s<-"Melanitta americana"
s<-"Glyptemis insculpta"
x<-as.data.table(resolve(s))
x <- gnr_resolve(s)
x <- as.data.table(name_backbone(name=s, rank='species'))
x<- name_lookup(query=s)

key<-name_backbone(s)$usageKey
x<-name_usage(key=key,data="synonyms")$data$canonicalName
x
x[x%in%ed$species]

as.data.table(name_backbone(s, rank='species'))
x<-as.data.table(name_lookup(query=s)$data)
x[,c("scientificName","species","canonicalName","accepted","rank")]
#grep(" x ",as.data.table(name_lookup(query=s)$data)$species)

w<-which(!ed$species%in%d$species)

l<-lapply(seq_along(w),function(i){
  print(i)
  x<-name_backbone(name=ed$species[w[i]], rank='species')

  cbind(ed$species[w[i]],x[,c("scientificName","species","synonym")])
})
x<-rbindlist(l)


ed$species[]


#######################################
## K, L and G functions ###############
#######################################

library(spatstat)

coords<-c("x","y")
win<-st_bbox(st_as_sf(obs[,..coords],coords=coords))
win<-owin(win[c(1,3)],win[c(2,4)])
pp<-as.ppp(unique(as.matrix(spobs[,..coords])),win)
reffort<-rasterize(as.matrix(obs[,..coords]),r[[1]],fun="length",background=0)
reffort<-aggregate(reffort,50,fun=sum,na.rm=TRUE)
plot(reffort)
effort<-as.im(list(x=xFromCol(reffort,1:ncol(reffort)),y=rev(yFromRow(reffort,1:nrow(reffort))),z=t(as.matrix(reffort,wide=TRUE)[nrow(reffort):1,])))
plot(effort,axes=TRUE)
effort<-list(effort=log(effort+0.1))
plot(effort$effort)
plot(pp,col="white",add=TRUE,cex=0.3,pch=16)

pp<-as.ppp(unique(as.matrix(obs[sample(1:nrow(obs),1000),..coords])),win)
fit<-ppm(pp, ~ effort, covariates = effort)
plot(predict(fit))
plot(pp,col="white",add=TRUE,cex=0.3,pch=16)
sims<-simulate(fit,drop=TRUE)
#fit<-ppm(sims, ~ effort, covariates = effort)
eem <- envelope(pp,Kinhom,funargs = list(lambda=fit),global=TRUE,nsim=19,fix.n=TRUE)
plot(eem)

robs<-rasterize(st_coordinates(st_as_sf(pp)[-1,]),reffort,fun="length",background=0)
observed<-values(robs)[,1]/res(robs)[1]^2
efforts<-log(values(reffort)[,1]+0.1)
plot(effectfun(fit,"effort"),ylim=c(0,max(observed)))
#robs<-rasterize(as.matrix(spobs[,..coords]),reffort,fun="length",background=0)
points(efforts,observed)


#######################################
#######################################
## Visual method ######################
#######################################
#######################################

library(FRutils)
library(magick)

inla.mesh2sp <- function(mesh,crs=NULL) {
#crs <- inla.CRS(inla.CRSargs(mesh$crs))
#isgeocentric <- identical(inla.as.list.CRS(crs)[["proj"]], "geocent")
#if (isgeocentric || (mesh$manifold == "S2")) {
#stop(paste0(
#"'sp' doesn't support storing polygons in geocentric #coordinates.\n",
#"Convert to a map projection with inla.spTransform() before
#calling inla.mesh2sp()."))
#}

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

grid<-function(im){
  inf<-image_info(im)
  w<-seq(0,inf$width,by=50)
  h<-seq(0,inf$height,by=50)
  for(i in seq_along(w)){
    im<-image_annotate(im,w[i],size=30,gravity="northwest",location=paste0("+",w[i],"+",inf$height*0.10),degrees=-90,weight=500)
  }
    for(i in seq_along(h)){
    im<-image_annotate(im,h[i],size=30,gravity="northwest",location=paste0("+",inf$width*0.01,"+",h[i]),degrees=0,weight=500)
  }
  im <- image_draw(im)
    abline(v = w, col = 'grey', lwd = '1', lty = 3)
    abline(h = h, col = 'grey', lwd = '1', lty = 3)
  dev.off()
  im<-image_border(im,"grey","2x2")
  im
}
#plot(grid(im))

### params
sp<-"Setophaga pinus"
crsr<-crs(r)
be<-unlist(ed[match(sp,ed$scientific_name),c("start","end")])
obs<-d[md>=be[1] & md<=be[2],]
obs<-unique(obs,by=c("recordedBy","species","cell"))
coords<-c("x","y")
obs<-st_as_sf(obs,coords=coords,crs=crsr)
spobs<-obs[obs$species==sp,]
# seagreen # orangered2 # snow2
cols<-list(
  all=unname(palette.colors()[6]),
  sp=unname(palette.colors()[7]),
  land="snow2"
)
cols$gradientN<-c("lightblue1","lightblue","cornflowerblue",cols$all,"darkblue","grey20")
cols$gradientPred<-c("grey90","lightgreen","forestgreen","black")
cols$gradientSDM<-c(cols$land,cols$sp,"darkred","darkred","grey10")

### 1 
png("/data/sdm_rbq/graphics/vabstract1.png",width=8,height=8,units="in",res=500,pointsize=11)
layout(1)
par(mar=c(0,0,0,0))
plot(st_geometry(na),col=cols$land,border="white")
plot(st_geometry(nalakes),col="white",border=adjustcolor("black",0.5),lwd=0.1,add=TRUE)
plot(st_geometry(obs),col=adjustcolor(cols$all,0.05),pch=16,cex=1,add=TRUE)
plot(st_geometry(spobs),col=adjustcolor(cols$sp,0.95),pch=16,cex=1.5,add=TRUE)
legend("topright",legend=c(sp,"All bird species"),pch=16,cex=2.1,col=rev(c(cols$all,cols$sp)),bty="n",inset=c(0,0.04),pt.cex=c(1.5,1),text.font=2)
plot(st_geometry(na),border=adjustcolor("black",0.5),lwd=1,add=TRUE)
dev.off()


### 2
ch<-rangeHull(spobs,breaks=100)
areas<-st_area(ch$hull)
png("/data/sdm_rbq/graphics/vabstract2_1.png",width=9,height=8,units="in",res=200,pointsize=11)
par(mar=c(4,4.5,1,1))
ylim<-c(0,max(c(ch$a,max(areas)))*1.05)
plot(ch$n,areas,ylim=ylim,yaxt="n",xaxt="n",xlab="",ylab="",col="black",pch=16,cex=1.5)
axis(1,tcl=-0.3,mgp=c(2,0.5,0),cex.axis=1.5,lwd=0,lwd.ticks=1)
axis(2,las=2,tcl=-0.3,mgp=c(2,0.5,0),cex.axis=1.5,lwd=0,lwd.ticks=1,at=pretty(areas),labels=pretty(areas)/10^6)
abline(h=ch$a,lty=3,lwd=5,col="black")
vs<-0:max(ch$n)
lines(vs,predict(ch$model,data.frame(X=vs)),lwd=4,col="black")
box(col=cols$land,lwd=5)
mtext(side=1,line=2.5,cex=2.5,font=2,text="Number of locations sampled")
mtext(side=2,line=2,cex=2.5,font=2,text=expression(bold("Convex Hull Area "%*%~10^6~km^2)))
text(0,ylim[2]*1.02,label="Asymptote",adj=c(0,1),cex=2.5,font=2)
dev.off()

png("/data/sdm_rbq/graphics/vabstract2_2.png",width=8,height=8,units="in",res=200,pointsize=11)
plot(st_geometry(na),col=cols$land,border=adjustcolor("black",0.5),lwd=1)
plot(st_geometry(nalakes),col="white",border=adjustcolor("black",0.5),lwd=0.1,add=TRUE)
plot(st_geometry(spobs),col=adjustcolor(cols$sp,0.95),pch=16,cex=1,add=TRUE)
plot(st_geometry(ch$hull[st_geometry_type(ch$hull)!="POINT"]),border=adjustcolor("black",0.15),lwd=2,add=TRUE)
dev.off()

im1<-image_read("/data/sdm_rbq/graphics/vabstract2_1.png")
im2<-image_trim(image_read("/data/sdm_rbq/graphics/vabstract2_2.png"))
im<-image_composite(im1, image_scale(im2,"x1000"),gravity="southeast",offset="+110+160")
image_write(im,"/data/sdm_rbq/graphics/vabstract2.png")



### 3
region<-st_buffer(concaveman(st_cast(na,"MULTIPOINT"),concavity=2),50)
region<-st_transform(region,crs(predictors))
set.seed(1234)
domain<-st_coordinates(st_sample(region,5000))
domain <- inla.nonconvex.hull(domain,convex = -0.015,resolution=75)
pedge<-0.02
edge<-min(c(diff(st_bbox(region)[c(1,3)])*pedge,diff(st_bbox(region)[c(2,4)])*pedge))
edge
Mesh <- inla.mesh.2d(loc.domain = NULL, #coordinates(occsp),
                     max.edge = c(edge,edge*2),
                     cutoff = edge,
                     offset = c(edge,edge*2),
                     boundary=domain,
                     crs = st_crs(na)
)
mesh<-inla.mesh2sp(Mesh,crs=CRS(prj))
plan(sequential)
weights <- ppWeight(sPoly = region, mesh = Mesh)
dmesh<-attributes(weights)$dmesh
dmesh$areas<-weights
dmesh$dmesh<-dmesh$id

png("/data/sdm_rbq/graphics/vabstract3_1.png",width=8,height=8,units="in",res=500,pointsize=11)
layout(matrix(1:,ncol=1,byrow=FALSE))
par(mar=c(0,0,0,0))
hublot<-st_buffer(st_transform(st_sfc(st_point(c(-74,44)),crs=4326),crs=st_crs(na)),dist=400)
xlim<-st_bbox(dmesh)[c(1,3)]
ylim<-st_bbox(dmesh)[c(2,4)]
plot(st_geometry(na),col=cols$land,border="black",xlim=xlim,ylim=ylim,lwd=0.1)
plot(st_geometry(nalakes),col="white",border="black",lwd=0.1,add=TRUE)
plot(mesh$triangles,border=adjustcolor("black",0.15),lty=3,add=TRUE)
plot(st_geometry(dmesh),border=adjustcolor("seagreen",0.25),add=TRUE)
plot(st_geometry(na),col=NA,border=adjustcolor("black",0.99),lwd=1,add=TRUE)
plot(st_geometry(hublot),border="forestgreen",lwd=5,add=TRUE)
legend("topright",legend=c("INLA mesh","Dual mesh"),lwd=c(1,4),lty=c(3,1),col=c(adjustcolor("black",0.85),adjustcolor("seagreen",0.35)),bty="n",cex=2.25,inset=c(0.0,0.01),seg.len=0.9)
dev.off()

png("/data/sdm_rbq/graphics/vabstract3_2.png",width=8,height=8,units="in",res=500,pointsize=11)
par(mar=c(0,0,0,0))
#xlim<-st_bbox(hublot)[c(1,3)]
#ylim<-st_bbox(hublot)[c(2,4)]
ring<-st_difference(st_buffer(st_geometry(hublot),dist=20),hublot)
plot(st_geometry(ring))
plot(st_geometry(na),col=cols$land,border=adjustcolor("black",0.25),add=TRUE)
plot(st_geometry(nalakes),col="white",border=adjustcolor("black",0.25),add=TRUE)
plot(st_geometry(na),col=NA,border=adjustcolor("black",0.99),lwd=5,add=TRUE)
plot(mesh$triangles,border=adjustcolor("black",0.25),lty=3,lwd=2,add=TRUE)
plot(st_geometry(dmesh),border=adjustcolor("seagreen",0.25),lwd=8,add=TRUE)
plot(st_difference(st_buffer(hublot,300),hublot),col="white",border=NA,add=TRUE)
plot(st_geometry(ring),col="red",border=NA,add=TRUE)
plot(st_geometry(hublot),border="forestgreen",lwd=15,add=TRUE)
dev.off()

im1<-image_read("/data/sdm_rbq/graphics/vabstract3_1.png")
im2<-image_read("/data/sdm_rbq/graphics/vabstract3_2.png")
im2<-image_fill(im2, "#FFFFFF00", point = "+1+1", fuzz = 0)
im2<-image_trim(im2)
im2<-image_fill(im2, "#FFFFFFFF", point = paste0("+1+",image_info(im2)$height/2), fuzz = 0)
im<-image_composite(im1, image_scale(im2,"x1900"),gravity="southwest",offset="+110+160")
image_write(im,"/data/sdm_rbq/graphics/vabstract3.png")


### 4
library(FRutils)
dmesh$nbobs<-lengths(st_intersects(dmesh,obs))
#pals<-function(i){adjustcolor(colo.scale(i,gradient),0.95)}
#plot(dmesh["nbobs"],lwd=0.1,pal=pals,logz=TRUE,nbreaks=50,key.pos=4,reset=FALSE,las=2)
#plot(st_geometry(na),col=cols$land,border="white",add=TRUE)
#plot(dmesh["nbobs"],border=NA,pal=pals,logz=TRUE,nbreaks=50,add=TRUE)

.image_scale2<-
function (z, col, breaks = NULL, key.pos, add.axis = TRUE, at = NULL,labels = NULL, 
    ..., axes = FALSE, key.length, logz = FALSE) 
{
    if (!is.null(breaks) && length(breaks) != (length(col) + 
        1)) 
        stop("must have one more break than colour")
    zlim = range(z, na.rm = TRUE)
    if (is.null(breaks)) 
        breaks = seq(zlim[1], zlim[2], length.out = length(col) + 
            1)
    if (is.character(key.length)) {
        kl = as.numeric(gsub(" cm", "", key.length))
        sz = if (key.pos %in% c(1, 3)) 
            dev.size("cm")[1]
        else dev.size("cm")[2]
        key.length = kl/sz
    }
    if (is.null(at)) {
        br = range(breaks)
        at = pretty(br)
        at = at[at > br[1] & at < br[2]]
    }
    kl_lim = function(r, kl) {
        m = mean(r)
        (r - m)/kl + m
    }
    if (key.pos %in% c(1, 3)) {
        ylim = c(0, 1)
        xlim = kl_lim(range(breaks), key.length)
        mar = c(0, ifelse(axes, 2.1, 1), 0, 1)
    }
    if (key.pos %in% c(2, 4)) {
        ylim = kl_lim(range(breaks), key.length)
        xlim = c(0, 1)
        mar = c(ifelse(axes, 2.1, 4), 0, 4, 0)
    }
    mar[key.pos] = 4.5
    par(mar = mar)
    poly = vector(mode = "list", length(col))
    for (i in seq(poly)) poly[[i]] = c(breaks[i], breaks[i + 
        1], breaks[i + 1], breaks[i])
    plot(1, 1, t = "n", ylim = ylim, xlim = xlim, axes = FALSE, 
        xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    offset = 0.01
    offs = switch(key.pos, c(0, 0, -offset, -offset), c(0, 0, 
        -offset, -offset), c(offset, offset, 0, 0), c(offset, 
        offset, 0, 0))
    for (i in seq_along(poly)) {
        if (key.pos %in% c(1, 3)) 
            polygon(poly[[i]], c(0, 0, 1, 1) + offs, col = col[i], 
                border = NA)
        if (key.pos %in% c(2, 4)) 
            polygon(c(0, 0, 1, 1) + offs, poly[[i]], col = col[i], 
                border = NA)
    }
    bx = c(breaks[1], rep(tail(breaks, 1), 2), breaks[1])
    #if (key.pos %in% c(1, 3)) 
    #    polygon(bx, c(0, 0, 1, 1) + offs, col = NA, border = cols$land,lwd=15.1)
    #if (key.pos %in% c(2, 4)) 
    #    polygon(c(0, 0, 1, 1) + offs, bx, col = NA, border = cols$land,lwd=15.1)
    #labels = if (logz) 
    #    parse(text = paste0("10^", at))
    #else if (inherits(breaks, c("POSIXt", "Date"))) 
    #    format(at)
    #else TRUE
    if (add.axis) 
        axis(key.pos, at = at, labels = labels, las=2, lwd=0, lwd.ticks=1,tcl=-0.2,hadj=0.2,cex.axis=3)
}

png("/data/sdm_rbq/graphics/vabstract4.png",width=9,height=8,units="in",res=500,pointsize=11)
layout(matrix(c(rep(1,91),rep(2,9)),ncol=100))
par(mar=c(0,0,0,0))
nbreaks<-100
legcolo<-exp(seq(log(1),log(max(dmesh$nbobs)),length.out=nbreaks))
vals<-c(dmesh$nbobs,legcolo)
vals2<-ifelse(vals==0,NA,log(vals))
colo<-colo.scale(vals2,cols$gradientN)
colo<-ifelse(vals==0,cols$land,colo)
plot(st_geometry(dmesh),col=colo[1:nrow(dmesh)],lwd=0.1,border=adjustcolor("black",0.5))
plot(st_geometry(dmesh),border=adjustcolor("black",0.15),add=TRUE)
plot(st_geometry(na),col=NA,border=adjustcolor("black",0.99),lwd=1,add=TRUE)
#plot(st_buffer(st_geometry(dmesh),dist=-20),col=colo[1:nrow(dmesh)],lwd=0.1,border=adjustcolor("black",0.5),add=TRUE)
#plot(st_geometry(na),col=NA,border=adjustcolor("black",0.25),lwd=1,add=TRUE)
#text(st_coordinates(st_centroid(dmesh)),label=dmesh$nbobs,cex=0.25)
at<-seq(log(1),log(max(dmesh$nbobs)),length.out=7)
labels<-round(exp(at))
at<-log(10^c(0:5))
labels<-parse(text = paste0("10^", 0:5))
.image_scale2(tail(log(legcolo),nbreaks),col=tail(colo,nbreaks),key.pos=4,key.length=0.8,logz=FALSE,at=at,labels=labels)
mtext(side=3,line=-4,adj=0.9,text="Number of\nobservations",cex=2.25,font=2)
dev.off()



### 5

# \overbrace{\sigma} = \underbrace{N} \times \overbrace{\left(\frac{r-r_{min}}{r_{max}-r_{min}} \times (S_{max}-1) + 1\right)}
# r = p/S
# p = Presence of target species (present = 1, absent = 0)
# S = Number of species observed

im<-image_read("/data/sdm_rbq/graphics/spspecificeff.png")
im<-image_border(im,"white","50x200")
im<-image_annotate(im,"Effort",size=45,gravity="northwest",location="+70+240",degrees=0,weight=700)
im<-image_annotate(im,"Number of observations\n      from target group",size=45,gravity="northwest",location="+115+465",degrees=0,weight=700)
im<-image_annotate(im,"Adjustement for species specific effort",size=45,gravity="northwest",location="+775+150",degrees=0,weight=700)
im<-image_trim(im)
im1<-image_border(im,"white","x50")

im<-image_read("/data/sdm_rbq/graphics/ratioeff.png")
im<-image_trim(im)
im2<-image_border(im,"white","x50")

im<-image_read("/data/sdm_rbq/graphics/peff.png")
im<-image_trim(im)
w<-image_info(im)$width+1500+40
im<-image_border(im,"white","1500")
im<-image_annotate(im,"Presence of target species (present = 1, absent = 0)",size=45,gravity="west",location=paste0("+",w,"-5"),degrees=0,weight=700)
im<-image_trim(im)
im3<-image_border(im,"white","x20")

im<-image_read("/data/sdm_rbq/graphics/Seff.png")
im<-image_trim(im)
w<-image_info(im)$width+1500+40
im<-image_border(im,"white","1500")
im<-image_annotate(im,"Number of species observed",size=45,gravity="west",location=paste0("+",w,"+8"),degrees=0,weight=700)
im<-image_trim(im)
im4<-image_border(im,"white","x20")

w1<-image_info(im1)$width
im<-image_append(c(im2,im3,im4),stack=TRUE)
w2<-image_info(im)$width
w<-(w1-w2)/2
im<-image_border(im,"white",paste(w))

im<-image_append(c(im1,im),stack=TRUE)
im<-image_border(im,"white","25x25")
im<-image_border(im,cols$land,"10x10")
im1<-image_border(im,"white","50x80")
#im1<-image_annotate(im,"Effort function",size=85,gravity="northwest",location="+80+0",degrees=0,weight=700)

#im<-image_read_pdf("/data/sdm_rbq/graphics/tableeff.pdf",density=600)
#im<-image_trim(im,fuzz=50)
#plot(im)

library(ggplot2)
library(gridExtra)
rescale<-function(x,to=c(1,2)){
  ((x-min(x))/(max(x)-min(x)))*(to[2]-to[1])+to[1]
}
effort<-function(n,N,S){
  p<-ifelse(n>0,1,0)
  r<-p/S
  #N*((((r-min(r))/(max(r)-min(r)))*(max(S)-1))+1)
  N*rescale(r,to=c(1,max(S)))
}
n<-c(0,0,0,0,0,0,10,1,1,5,5,5)
N<-c(1,1,1,10,10,10,10,10,10,100,100,100)
S<-c(1,1,1,1,5,10,1,5,10,10,50,96)
Smax<-max(S)
n_N<-formatC(n/N,digits=3,format="f")
n_N<-paste(n,N,sep="/")
n_N<-n/N
eg<-data.frame(n,N,S,n_N,Smax)
eg$sigma<-round(effort(n,N,S),1)
fills<-rep(c("grey90","grey80","grey90","grey80"),each=3)
tt <- ttheme_default(
  colhead=list(fg_params=list(parse=TRUE,fontface=4),bg_params=list(fill=cols$all)),
  core=list(bg_params=list(fill=fills, col=NA))
)
headers<-c('n','N','S','n / N','Smax','sigma')
g<-tableGrob(eg,rows=NULL,cols=headers,widths=unit(rep(1,ncol(eg)),c("in")),theme=tt)
#grid.arrange(g,widths=20)
ggsave(plot=g, filename="/data/sdm_rbq/graphics/tableeff.png")
im2<-image_read("/data/sdm_rbq/graphics/tableeff.png")
im2<-image_trim(im2)
im2<-image_border(im2,"white","10x10")
image_write(im2,"/data/sdm_rbq/graphics/tableeff.png")
im2<-image_read("/data/sdm_rbq/graphics/tableeff.png")



im2<-image_read("/data/sdm_rbq/graphics/tableeff.png")
im2<-image_trim(im2)
w1<-image_info(im1)$width
w2<-image_info(im2)$width
w<-(w1-w2)/2
im2<-image_border(im2,"white",w)

im<-image_append(c(im1,im2),stack=TRUE)

png("/data/sdm_rbq/graphics/vabstract5.png",width=8,height=8,units="in",res=500,pointsize=11)
par(mar=c(0,0,0,0),bg="white")
plot(im)
dev.off()


 #rescale(c(0/20,1/20,20/20),to=c(1,20))





### 6
op<-rast("/data/predictors_sdm/predictors.tif")
xv<-"forested"
rop<-op[xv]
rop<-extend(rop,extend(ext(dmesh),res(rop)))
xy<-xyFromCell(rop,1:ncell(rop))
xy<-cbind(xy,values(rop))
xynotna<-xy |> as.data.table() |> na.omit() |> as.matrix()
nn<-knnx.index(xynotna[,1:2],xy[,1:2],k=1)
rop<-setValues(rop,xynotna[nn,-(1:2)])
e<-exact_extract(rop[xv], 
                        dmesh, 
                        fun = function(values, coverage_fraction){
                        colSums(as.matrix(values) * coverage_fraction,
                                na.rm = TRUE) / sum(coverage_fraction)
                        },
                        force_df = FALSE,
                        progress = TRUE)
dmesh[[xv]]<-e

png("/data/sdm_rbq/graphics/vabstract6.png",width=9,height=8,units="in",res=500,pointsize=11)
layout(matrix(c(rep(1,91),rep(2,9)),ncol=100))
par(mar=c(0,0,0,0))
plot(st_geometry(dmesh),col=colo.scale(dmesh[[xv]],cols$gradientPred),mar=c(0,0,0,0),border="black",lwd=0.1)
plot(st_geometry(dmesh),border=adjustcolor("black",0.15),add=TRUE)
plot(st_geometry(na),col=NA,border=adjustcolor("black",0.99),lwd=1,add=TRUE)
vals<-seq(min(dmesh[[xv]],na.rm=TRUE),max(dmesh[[xv]],na.rm=TRUE),length.out=200)
.image_scale2(vals,col=colo.scale(vals,cols$gradientPred),key.pos=4,key.length=0.8)
mtext(side=3,line=-4,adj=0.9,text="Proportion\nof forests",cex=2.25,font=2)
dev.off()


### 7
# https://latex2png.com/
im<-image_read("/data/sdm_rbq/graphics/lgcpeff.png")
#im<-image_border(im,"red","5x5")
im<-image_border(im,"white","50x200")
im<-image_annotate(im,"Intensity",size=45,gravity="northwest",location="+40+240",degrees=0,weight=700)
im<-image_annotate(im,"Effort",size=45,gravity="northwest",location="+332+495",degrees=0,weight=700)
im<-image_annotate(im,"Intercept + predictors",size=45,gravity="northwest",location="+460+150",degrees=0,weight=700)
im<-image_annotate(im,"Spatial field",size=45,gravity="northwest",location="+895+405",degrees=0,weight=700)
im<-image_trim(im)
im<-image_border(im,"white","50x50")
im<-image_border(im,cols$land,"10x10")
im<-image_scale(im,image_info(im2)$width)
b<-(image_info(im2)$width-image_info(im)$width)/2
im<-image_border(im,"white",b)
im<-image_border(im,"white","25x150")
im<-image_annotate(im,"Log Gaussian Cox process adjusted for effort",size=150,gravity="northwest",location="+80+0",degrees=0,weight=700)
im<-image_border(im,"white","50x100")
#image_write(im,"/data/sdm_rbq/graphics/vabstract7.png")
png("/data/sdm_rbq/graphics/vabstract7.png",width=8,height=8,units="in",res=500,pointsize=11)
par(mar=c(0,0,0,0),bg="white")
plot(im)
dev.off()


### 8
png("/data/sdm_rbq/graphics/vabstract8_1.png",width=9,height=8,units="in",res=500,pointsize=11)
sdm<-rast(paste0("/data/sdm_rbq/rasters/",gsub(" ","_",sp),"_birds.tif"))$mean
layout(matrix(c(rep(1,91),rep(2,9)),ncol=100))
par(mar=c(0,0,0,0))
plot(sdm,col=colo.scale(1:200,cols$gradientSDM),axes=FALSE,mar=c(0,0,0,0),legend=FALSE)
plot(st_geometry(na),col=NA,border=adjustcolor("black",0.5),lwd=1,add=TRUE)
plot(st_geometry(nalakes),col="white",border=adjustcolor("black",0.5),lwd=0.1,add=TRUE)
vals<-seq(min(values(sdm),na.rm=TRUE),max(values(sdm),na.rm=TRUE),length.out=200)
.image_scale2(vals,col=colo.scale(vals,cols$gradientSDM),key.pos=4,key.length=0.8)
mtext(side=4,line=-4,adj=0.15,text="Intensity",cex=2.25,font=2)
dev.off()

png("/data/sdm_rbq/graphics/vabstract8_2.png",width=9,height=8,units="in",res=500,pointsize=11)
sdm<-rast(paste0("/data/sdm_rbq/rasters/",gsub(" ","_",sp),"_birds.tif"))$linksd
layout(matrix(c(rep(1,91),rep(2,9)),ncol=100))
par(mar=c(0,0,0,0))
plot(sdm,col=colo.scale(1:200,cols$gradientSDM),axes=FALSE,mar=c(0,0,0,0),legend=FALSE)
plot(st_geometry(na),col=NA,border=adjustcolor("black",0.5),lwd=1,add=TRUE)
plot(st_geometry(nalakes),col="white",border=adjustcolor("black",0.5),lwd=0.1,add=TRUE)
vals<-seq(min(values(sdm),na.rm=TRUE),max(values(sdm),na.rm=TRUE),length.out=200)
.image_scale2(vals,col=colo.scale(vals,cols$gradientSDM),key.pos=4,key.length=0.8)
mtext(side=3,line=-4,adj=0.85,text="SD (link scale)",cex=3.5,font=2)
dev.off()

im<-image_read("/data/sdm_rbq/graphics/vabstract8_1.png")
im2<-image_trim(image_read("/data/sdm_rbq/graphics/vabstract8_2.png"))
im2<-image_fill(im2, "#FFFFFF00", point = "+18+18", fuzz = 0)
im<-image_composite(im, image_scale(im2,"x1000"),gravity="northeast",offset="+500+120")
image_write(im,"/data/sdm_rbq/graphics/vabstract8.png")


############################################
############################################
############################################

framecol<-"lemonchiffon2"

titles<-c("A. Multispecies occurrences\nwith target species","B. Assess coverage using\nstabilisation of convex hull","C. Build mesh\ncovering study area","D. Summarize N for each\ncell of the dual mesh","E. Adjust N to reflect\nspecies-specific effort","F. Extend and summarize\npredictors to dual mesh area","G. Combine model\ncomponents and run INLA","H. Map intensity along\nwith uncertainty")
images<-list.files("/data/sdm_rbq/graphics/",pattern="vabstract",full.names=TRUE)
images<-images[-grep("_",basename(images))]
images<-images[rep(1:length(images),length.out=length(titles))]

ims<-do.call("c",lapply(seq_along(images),function(x){
  im<-image_read(images[x])
  im<-image_scale(im,"x700")
  im<-image_border(im,framecol,"5x5")
  im<-image_border(im,framecol,"25x125")
  im<-image_annotate(im,titles[x],size=50,gravity="northwest",location="+40+30",degrees=0,weight=900)
  im<-image_trim(im)
  im<-image_border(im,framecol,"25x25")
  #im<-image_border(im,"grey90","5x5")
  #im<-image_fill(im, framecol, point = "+5+5", fuzz = 0)
  #im<-image_fill(im, framecol, point = "+1+1", fuzz = 0)
  #im<-image_border(im, framecol,"10x10")
  im<-roundCorners(im,rounding=0.05,corners=1:4)
  #im<-image_trim(im)
  #im<-image_fill(im, framecol, point = "+21+21", fuzz = 0)
  im<-image_border(im,"white","25x25",operator="atop")
  im
}))

im1<-image_append(ims[1:4], stack = FALSE)
#im1<-image_border(ims[5:8],"white","x10")
im2<-image_append(ims[5:8], stack = FALSE)
im<-image_append(c(im1,im2), stack = TRUE)
image_write(im,"/data/sdm_rbq/graphics/visual_abstract.png")



#print(image_info(im))
#im <- image_draw(im)
#arrows(x0=1300, y0=600, x1 = 1400, y1 = 600, length = 0.3, angle = 30, lwd=40,code=2,col="seagreen")
#dev.off()
#image_write(im,"/data/sdm_rbq/graphics/test.png")


p<-as.integer(n>0)
N*scales::rescale(n/N,to=c(1,max(S)))
N*scales::rescale(p/S,to=c(1,max(S)))


library(interp)
IL <- interp(xi,yi,fi,xo,yo,method="linear")

k<-1:nrow(dmesh)
#mpp$summary.fitted.values[k,"mean"]

pr<-aggregate(r[[1]],100)
xy<-xyFromCell(pr,1:ncell(pr))
o<-st_intersects(st_as_sf(data.frame(xy),coords=c("x","y"),crs=st_crs(na)),na)
xy<-xy[lengths(o)>0,]
test <- interp(x=Mesh$loc[,1],y=Mesh$loc[,2],z=mpp$summary.fitted.values[k,"mean"],xo=xy[,1],yo=xy[,2],method="linear",output="points")


test <- interp(x=Mesh$loc[,1],y=Mesh$loc[,2],z=mpp$summary.fitted.values[k,"mean"],xo=xy[,1],yo=xy[,2],method="linear",output="grid")


#####################################################
#####################################################
### Species maps ####################################
#####################################################
#####################################################

library(magick)
library(berryFunctions)
library(ebirdst)
library(terra)
library(sf)
library(FRutils)
library(ebirdst)
library(data.table)
library(rgbif)

Sys.setlocale("LC_TIME", "English")

pics<-list(
  Setophaga_pinus="https://static.inaturalist.org/photos/115485881/large.jpg",
  Setophaga_petechia="https://inaturalist-open-data.s3.amazonaws.com/photos/191200849/large.jpg",
  Setophaga_ruticilla="https://inaturalist-open-data.s3.amazonaws.com/photos/31539046/large.jpg",
  Setophaga_americana="https://inaturalist-open-data.s3.amazonaws.com/photos/135474901/large.jpg",
  Setophaga_citrina="https://inaturalist-open-data.s3.amazonaws.com/photos/19448070/large.jpg",
  Setophaga_virens="https://inaturalist-open-data.s3.amazonaws.com/photos/199207488/large.jpeg"
)

path<-"C:/Users/God/Downloads"
colmean<-c("#DDDEE0","#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00","darkorange")
colmean<-c("#CCCCCC","#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00","darkorange")

ed<-as.data.table(ebirdst_runs)
ed<-ed[,.(common_name,scientific_name,species_code,resident,breeding_start,breeding_end)]
startend<-c("breeding_start","breeding_end")
ed[ ,(startend):=lapply(.SD,as.character),.SDcols=startend]
ed[,breeding_start:=ifelse(is.na(breeding_start),"2020-01-01",breeding_start)]
ed[,breeding_end:=ifelse(is.na(breeding_end),"2020-12-31",breeding_end)]
ed[,start:=substr(breeding_start,6,10)]
ed[,end:=substr(breeding_end,6,10)]
ed[,species:=scientific_name]

sdm<-rast(paste0("C:/Users/God/Downloads/","Setophaga_americana","_birds.tif"))
na<-st_read(list.files("C:/Users/God/Documents/predictors_sdm",pattern="na.shp",full=TRUE))
na<-st_transform(na,st_crs(sdm))
coast<-st_union(na)
nalakes<-st_read("C:/Users/God/Documents/predictors_sdm/nalakes.gpkg")
nalakes<-st_transform(nalakes,st_crs(na))
areas<-as.numeric(st_area(nalakes))
ke<-paste(c("Manicouagan","-Jean"),collapse="|")
nalakes<-nalakes[unique(c(rev(order(areas))[1:50],grep(ke,nalakes$name_fr))),]

x<-do.call("rbind",lapply(names(pics),function(i){
  o<-occ_search(scientificName=gsub("_"," ",i),hasCoordinate=TRUE,limit=1000,month="6,7")
  keep<-c("species","decimalLongitude","decimalLatitude")
  x<-o$data[,keep]
  x<-st_as_sf(x,coords=keep[2:3],crs=4326)
  x<-st_transform(x,st_crs(na))
  x$sp<-gsub(" ","_",x$species)
  x
}))
x<-x[na,]

.image_scale2<-
  function (z, col, breaks = NULL, key.pos, add.axis = TRUE, at = NULL,labels = NULL, 
            ..., axes = FALSE, key.length, logz = FALSE) 
  {
    if (!is.null(breaks) && length(breaks) != (length(col) + 
                                               1)) 
      stop("must have one more break than colour")
    zlim = range(z, na.rm = TRUE)
    if (is.null(breaks)) 
      breaks = seq(zlim[1], zlim[2], length.out = length(col) + 
                     1)
    if (is.character(key.length)) {
      kl = as.numeric(gsub(" cm", "", key.length))
      sz = if (key.pos %in% c(1, 3)) 
        dev.size("cm")[1]
      else dev.size("cm")[2]
      key.length = kl/sz
    }
    if (is.null(at)) {
      br = range(breaks)
      at = pretty(br)
      at = at[at > br[1] & at < br[2]]
    }
    kl_lim = function(r, kl) {
      m = mean(r)
      (r - m)/kl + m
    }
    if (key.pos %in% c(1, 3)) {
      ylim = c(0, 1)
      xlim = kl_lim(range(breaks), key.length)
      mar = c(0, ifelse(axes, 2.1, 1), 0, 1)
    }
    if (key.pos %in% c(2, 4)) {
      ylim = kl_lim(range(breaks), key.length)
      xlim = c(0, 1)
      mar = c(ifelse(axes, 2.1, 4), 0, 4, 0)
    }
    mar[key.pos] = 4.5
    par(mar = mar)
    poly = vector(mode = "list", length(col))
    for (i in seq(poly)) poly[[i]] = c(breaks[i], breaks[i + 
                                                           1], breaks[i + 1], breaks[i])
    plot(1, 1, t = "n", ylim = ylim, xlim = xlim, axes = FALSE, 
         xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    offset = 0.01
    offs = switch(key.pos, c(0, 0, -offset, -offset), c(0, 0, 
                                                        -offset, -offset), c(offset, offset, 0, 0), c(offset, 
                                                                                                      offset, 0, 0))
    for (i in seq_along(poly)) {
      if (key.pos %in% c(1, 3)) 
        polygon(poly[[i]], c(0, 0, 1, 1) + offs, col = col[i], 
                border = NA)
      if (key.pos %in% c(2, 4)) 
        polygon(c(0, 0, 1, 1) + offs, poly[[i]], col = col[i], 
                border = NA)
    }
    bx = c(breaks[1], rep(tail(breaks, 1), 2), breaks[1])
    #if (key.pos %in% c(1, 3)) 
    #    polygon(bx, c(0, 0, 1, 1) + offs, col = NA, border = cols$land,lwd=15.1)
    #if (key.pos %in% c(2, 4)) 
    #    polygon(c(0, 0, 1, 1) + offs, bx, col = NA, border = cols$land,lwd=15.1)
    #labels = if (logz) 
    #    parse(text = paste0("10^", at))
    #else if (inherits(breaks, c("POSIXt", "Date"))) 
    #    format(at)
    #else TRUE
    if (add.axis) 
      axis(key.pos, at = at, labels = labels, las=2, lwd=0, lwd.ticks=1,tcl=-0.2,hadj=0.2,cex.axis=3)
  }


roundCorners<-function(im,rounding=0.05,corners=1:4){
  #im<-image_read("C:/Users/God/Downloads/Setophaga_petechia_map.png")
  #im<-image_scale(im,"300")
  w<-image_info(im)$width
  h<-image_info(im)$height
  mask <- image_graph(width = w, height = h, res = 96)
  par(mar=c(0,0,0,0))
  plot(0:1,0:1,type="n",xaxt="n",yaxt="n",bty="n",bg="magenta")
  dev.off()
  mask <- image_draw(mask)
  roundedRect(0,h,w,0,rounding=rounding,col="black",corners=corners)
  dev.off()
  res<-image_composite(mask, im, "plus")
  res
}
#roundCorners(image_scale(maps2,"x400"))


roundIm<-function(url,file,path,open=FALSE){
  mask <- image_graph(width = 400, height = 400, res = 96)
  par(mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",bty="n")
  dev.off()
  mask <- image_draw(mask)
  symbols(200,200,circles=198,bg=1,inches=FALSE,add=TRUE)
  dev.off()
  im <- image_read(url)
  #im<-image_scale(im,"200")
  wh<-c(image_info(im)$width,image_info(im)$height)
  mi<- min(wh)
  wm<-which.max(wh)
  if(wh[1]!=wh[2]){
    if(wm==1){
      geom<-paste0(mi, "x", mi, "+",abs(diff(wh))/2,"+","0")
    }else{
      geom<-paste0(mi, "x", mi, "+","0","+",abs(diff(wh))/2)
    }
  }
  im<-image_crop(im, geometry=geom,repage=TRUE)
  mask<-image_scale(mask, as.character(image_info(im)$height))
  organism<-image_composite(mask, im, "plus") 
  organism<-image_scale(organism,"500")
  organim<-image_trim(organism)
  h<-image_info(organism)$height
  organism<-image_border(organism,"#FFFFFF","50x50")
  organism<-image_fill(organism,"none",point="+5+5",fuzz=5)
  organism<-image_draw(organism)
  symbols(image_info(organism)$height/2,image_info(organism)$width/2,circles=h/2,bg="transparent",fg="grey99",inches=FALSE,lwd=20,add=TRUE)
  dev.off()
  organim<-image_trim(organism,fuzz=50)
  image_write(organism,file.path(path,paste0(file,".png")))
}

lapply(seq_along(pics),function(i){
  roundIm(url=pics[[i]],file=names(pics)[i],path=path,open=FALSE)
})



#im1<-image_read("C:/Users/God/Downloads/visual_abstract.png")
#im2<-image_read("C:/Users/God/Downloads/organism.png")
#im2<-image_scale(im2,"x225")
#im<-image_composite(im1, im2, gravity="northwest",offset="+20+350")
#image_write(im,"C:/Users/God/Downloads/compo.png")


lapply(names(pics),function(i){
  sdm<-rast(paste0("C:/Users/God/Downloads/",i,"_birds.tif"))
  png(file.path(path,paste0(i,"_map.png")),width=4,height=4,res=400,units="in")
  plot(sdm$mean,col=colo.scale(1:200,colmean),mar=c(0,6,0,0),bty="n",axes=FALSE,legend=FALSE)
  plot(st_geometry(na),border="white",lwd=0.1,add=TRUE)
  plot(st_geometry(nalakes),col="white",border=adjustcolor("black",0.15),lwd=0.1,add=TRUE)
  plot(st_geometry(coast),lwd=0.1,border=gray(0,0.45),add=TRUE)
  plot(st_geometry(x[x$sp==i,]),pch=16,lwd=0.25,cex=0.15,col="grey1",add=TRUE)
  dev.off()
  png(file.path(path,paste0(i,"_sd.png")),width=4,height=4,res=400,units="in")
  plot(sdm$linksd,col=colo.scale(1:200,colmean),mar=c(0,6,0,0),bty="n",axes=FALSE,legend=FALSE)
  plot(st_geometry(na),border="white",lwd=0.1,add=TRUE)
  plot(st_geometry(nalakes),col="white",border=adjustcolor("black",0.15),lwd=0.1,add=TRUE)
  plot(st_geometry(coast),lwd=0.1,border=gray(0,0.45),add=TRUE)
  plot(st_geometry(x[x$sp==i,]),pch=16,lwd=0.25,cex=0.15,col="grey1",add=TRUE)
  dev.off()
  map<-image_read(file.path(path,paste0(i,"_map.png")))
  map<-image_trim(map)
  b<-200
  map<-image_border(map,"white",b)
  map<-image_crop(map,paste0(image_info(map)$width-b,"x100000"))
  lsd<-image_read(file.path(path,paste0(i,"_sd.png")))
  lsd<-image_trim(lsd)
  lsd<-image_scale(lsd,"x400")
  lsd<-image_fill(lsd, "#FFFFFF00", point = "+5+5", fuzz = 0)
  map<-image_composite(map, lsd, gravity="west",offset="+0-60")
  im<-image_border(map,"white","20x20")
  hi<-image_trim(image_read(file.path(path,"overlap.png")))
  hi<-image_scale(hi,"x350")
  hi<-image_fill(hi,"#FFFFFF00", point = "+1+1", fuzz = 1)
  im<-image_composite(im, hi, gravity="southwest",offset="+15+15")
  im<-image_border(im,"white","20x20")
  sp<-image_read(file.path(path,paste0(i,".png")))
  sp<-image_scale(sp,"x350")
  #im<-image_border(im,"white","10x10")
  im<-image_composite(im, sp, gravity="northeast",offset="+1+10")
  im<-image_fill(im,"lemonchiffon2", point = "+5+5", fuzz = 1) # palegoldenrod or lavender
  text<-gsub("_"," ",i)
  im<-image_annotate(im,"SD",size=50,color="grey40",strokecolor="white",gravity="west",location="+40-100",degrees=0,weight=700)
  m<-match(text,ed$scientific_name)
  dates<-paste(format(as.Date(c(ed$breeding_start[m],ed$breeding_end[m])),"%b %d"),collapse=" - ")
  im<-image_annotate(im,paste(" ",paste(text,dates,sep=" / ")," "),size=60,color="grey1",strokecolor="white",gravity="northwest",location="+0+0",degrees=0,weight=700,boxcolor=adjustcolor("white",0.5))
  im<-image_annotate(im,"re = 0.94",size=60,color="grey40",strokecolor="white",gravity="northwest",location="+15+95",degrees=0,weight=700,boxcolor="none")
  im<-image_annotate(im,"n = 213",size=60,color="grey40",strokecolor="white",gravity="northwest",location="+15+163",degrees=0,weight=700,boxcolor="none")
  
  ### scale
  #vals<-seq(min(values(sdm$mean),na.rm=TRUE),max(values(sdm$mean),na.rm=TRUE),length.out=200)
  vals<-seq(min(values(sdm$mean),na.rm=TRUE),max(values(sdm$mean),na.rm=TRUE),length.out=200)
  vals<-vals/max(vals)
  png(paste0("C:/Users/God/Downloads/",i,"_scale.png"),width=1.5,height=20,units="in",res=400)
  par(bg="lemonchiffon2")
  .image_scale2(vals,colo.scale(vals,colmean),key.length=0.9,key.pos=4,add=TRUE)
  dev.off()
  sc<-image_read(paste0("C:/Users/God/Downloads/",i,"_scale.png"))
  sc<-image_trim(sc)
  sc<-image_scale(sc,"x500")
  #image_write(sc,paste0("C:/Users/God/Downloads/",i,"_scale.png"))
  #file.show(paste0("C:/Users/God/Downloads/",i,"_scale.png"))
  im<-image_composite(im,sc,gravity="southeast",offset="+30+50")
  
  im<-roundCorners(im,corners=c(1,3:4))
  im<-image_border(im,"white","20x20")
  image_write(im,file.path(path,paste0(i,"_map.png")))
  #file.show(file.path(path,paste0(i,"_maps.png")))
})


ims<-do.call("c",lapply(names(pics),function(i){
  image_read(file.path(path,paste0(i,"_map.png")))
}))
im1<-image_append(ims[1:3],stack=TRUE)
im2<-image_append(ims[4:6],stack=TRUE)
maps<-image_append(c(im1,im2),stack=FALSE)
maps<-image_border(maps,"white","20x20")
image_write(maps,file.path(path,"maps.png"))
file.show(file.path(path,"maps.png"))


#maps2<-image_fill(maps, "red", point = "+5+5", fuzz = 1)
#maps2<-image_fill(maps2, adjustcolor("khaki",0.95), point = "+80+80", fuzz = 1)
#maps2<-image_fill(maps, "antiquewhite", point = "+5+5", fuzz = 1)
#image_write(maps2,file.path(path,"maps2.png"))
#file.show(file.path(path,"maps2.png"))

#vals<-seq(min(values(sdm$mean),na.rm=TRUE),max(values(sdm$mean),na.rm=TRUE),length.out=200)
#png(paste0("C:/Users/God/Downloads/",i,"_scale.png"),width=2,height=30,units="in",res=400)
#par(bg="lemonchiffon2")
#.image_scale2(vals,colo.scale(vals,colmean),key.length=0.9,key.pos=4,add=TRUE)
#dev.off()
#sc<-image_read(paste0("C:/Users/God/Downloads/",i,"_scale.png"))
#sc<-image_trim(sc)
#image_write(sc,paste0("C:/Users/God/Downloads/",i,"_scale.png"))
#file.show(paste0("C:/Users/God/Downloads/",i,"_scale.png"))



#########################################################
#########################################################
#########################################################
### Side by side with eBird #############################
#########################################################
#########################################################
#########################################################

#list.files(tools::R_user_dir("ebirdst"),recursive=TRUE)
#options(vsc.dev.args = list(width = 1500, height = 800))
colmean<-c("#CCCCCC","#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00","darkorange")
colmean<-colo.scale(1:200,colmean)
coast<-st_union(na)

naplot<-function(){
  plot(st_geometry(na),border="white",lwd=0.25,add=TRUE)
  plot(st_geometry(nalakes),col="white",border=adjustcolor("black",0.15),lwd=0.25,add=TRUE)
  plot(st_geometry(coast),lwd=0.25,border=gray(0,0.45),add=TRUE)
}

df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
spc<-df$species[which(df$I>=0.75 & df$I<=0.95 & df$reach>=0.85 & df$n>=100)]
spc

lsp<-sort(unique(df$species))

#for(i in lsp[1:10]){
cl<-makeCluster(5)
registerDoParallel(cl)
foreach(i=lsp,.packages=c("sf","terra","magick")) %dopar% {
  print(i)
  sp<-i
  sdm<-rast(file.path("/data/sdm_rbq/rasters",paste0(gsub(" ","_",sp),"_birds.tif")))[[1]]
  ebirdpath<-file.path("/data/predictors_sdm/expert_maps/eBird/abundance",paste0(gsub(" ","_",sp),"_ebird2.tif"))
  if(file.exists(ebirdpath)){  
    ebird<-rast(ebirdpath)
  }else{
    return(NULL) #next
  }

  fname<-paste0(gsub(" ","_",sp),"_comp")

  png(file.path("/data/sdm_rbq/temp",paste0(fname,"1.png")),width=5,height=5,res=400,units="in")
  plot((sdm)^(1/1),col=colmean,mar=c(0,0,0,0),legend=FALSE,axes=FALSE)
  naplot()
  dev.off()
  png(file.path("/data/sdm_rbq/temp",paste0(fname,"2.png")),width=5,height=5,res=400,units="in")
  plot((ebird)^(1/2),col=colmean,mar=c(0,0,0,0),legend=FALSE,axes=FALSE)
  naplot()
  dev.off()

  im1<-image_trim(image_read(file.path("/data/sdm_rbq/temp",paste0(fname,"1.png"))))
  im1<-image_border(im1,"white","2000")
  im2<-image_trim(image_read(file.path("/data/sdm_rbq/temp",paste0(fname,"2.png"))))
  im2<-image_fill(im2,"#FFFFFF00","+1+1")
  im<-image_composite(im1,im2,gravity="center",offset="+1600+0")
  im<-image_trim(im)
  im<-image_border(im,"white","50x50")
  im<-image_annotate(im,sp,size=120,color="black",gravity="northwest",location="+10+10",boxcolor="#FFFFFFAA",degrees=0,weight=700)
  m<-match(sp,df$species)
  im<-image_annotate(im,paste0("n = ",df$n[m]),size=80,color="black",gravity="west",location="+10-100",degrees=0,weight=700)
  im<-image_annotate(im,paste0("reach = ",round(df$reach,2)[m]),size=80,color="black",gravity="west",location="+10+0",degrees=0,weight=700)
  im<-image_annotate(im,paste0("cor = ",round(df$correlation[m],2)),size=80,color="black",gravity="west",location="+10+100",degrees=0,weight=700)
  im<-image_annotate(im,paste0("I = ",round(df$I[m],2)),size=80,color="black",gravity="west",location="+10+200",degrees=0,weight=700)
  im<-image_annotate(im,"mapSpecies",size=100,color="#00000044",gravity="northwest",location="+750+300",boxcolor="#FFFFFFAA",degrees=0,weight=700)
  im<-image_annotate(im,"eBird",size=100,color="#00000044",gravity="northwest",location="+2550+300",boxcolor="#FFFFFFAA",degrees=0,weight=700)
  im<-image_scale(im,"x500")
  image_write(im,file.path("/data/sdm_rbq/comparison",paste0(gsub(" ","_",sp),"_comparison.png")))
}
stopCluster(cl)


##########################################################
##########################################################
### Stacked uncertainty ##################################
##########################################################
##########################################################

lsdms<-list.files("/data/sdm_rbq/rasters",pattern="_birds.tif",full=TRUE)
lsdms<-lsdms[rev(order(file.info(lsdms)$mtime))]
df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
sp<-gsub(" ","_",df$species[which(df$reach>=0.00)])
lsdms<-lsdms[gsub("_birds.tif","",basename(lsdms))%in%sp]
sdms<-rast(lapply(lsdms,rast,lyrs="linksd"))
gsd<-mean(sdms)
#global(sdms,"max",na.rm=TRUE)

#plot(gsd/global(gsd,"min",na.rm=TRUE)[1,1],col=colmean,mar=c(0,0,0,0))
png("/data/sdm_rbq/graphics/stacked_uncertainty.png",width=6,height=5,res=400,units="in")
plot(gsd,col=c(rev(magma(100)),gray((0:200)/200)),mar=c(0,0,0,4),bty="n",axes=FALSE,legend=TRUE)
naplot()
dev.off()