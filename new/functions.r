
#########################################################################
### Some useful functions ###############################################
#########################################################################

### Some command line tools
# nohup Rscript inat_sdm.R > output.out 2>&1 &
# cat output.out
# grep -n 10 output.out
# tail -n 10 output.out ; grep -n -i checkpoint output.out

# find /tmp -mmin +60 -user rouf1703 -delete
# find /tmp -mmin +60 -user rouf1703 -type d -empty -exec rmdir {} \;
# df -h
# ls -laht | sort -k4 | grep rouf1703

##########################################################################
### Prints time at different checkpoints to monitor progress of script

checkpoint <- function(msg = "", tz = "Indian/Mauritius") {
  cat(paste(
    "checkpoint",
    format(Sys.time(), tz = tz, usetz = FALSE),
    msg,
    sep = " - "
  ), "\n")
}

#########################################################################
### Build new data with two way interactions and polynomials

newdata2 <- function(x,
                     n = 10,
                     n2 = 3,
                     fun = mean) {
  vn <-
    unique(unlist(strsplit(gsub(
      "[[:digit:]]", "", names(x)
    ), "\\."))) # all unique variable names
  vs <-
    unique(unlist(strsplit(names(x), "\\."))) # all vars or components of interactions
  vi <- names(x)[grep("\\.", names(x))] # all interactions
  fix <-
    as.data.frame(as.list(apply(x[, vn, drop = FALSE], 2, FUN = fun))) # mean values of each var
  # build variables set (corresponds to all graphs)
  temp <- strsplit(names(x), "\\.")
  len <- sapply(temp, length)
  if (any(len > 1)) {
    ints <- temp[len > 1]
    ints <- lapply(ints, function(i) {
      gsub("[[:digit:]]", "", i)
    })
    ints <- ints[!duplicated(sapply(ints, paste0, collapse = ""))]
    varset <- c(setdiff(vn, unlist(ints)), ints)
  } else{
    varset <- vn
  }
  # loop through set
  ans <- lapply(varset, function(vars) {
    if (length(vars) == 1) {
      nd <-
        data.frame(seq(min(x[, vars[1]], na.rm = TRUE), max(x[, vars[1]], na.rm =
                                                              TRUE), length = n))
      names(nd) <- vars[1]
    } else{
      nd <- expand.grid(seq(min(x[, vars[1]], na.rm = TRUE), max(x[, vars[1]], na.rm =
                                                                   TRUE), length = n),
                        seq(min(x[, vars[2]], na.rm = TRUE), max(x[, vars[2]], na.rm = TRUE), length =
                              n2))
      names(nd) <- vars
    }
    # add variables not implied in interactions
    notvars <- setdiff(vn, vars)
    if (length(notvars)) {
      nd <- cbind(nd, fix[rep(1, nrow(nd)), notvars])
      names(nd) <- c(vars, notvars)
    }
    # add polynomials
    poly <-
      apply(expand.grid(vn, as.character(2:10)), 1, paste, collapse = "") # possibilities in polynomials (check up to 10)
    elev <- vs[vs %in% poly]
    if (length(elev)) {
      l <- lapply(elev, function(i) {
        la <- gsub("[[:digit:]]", "", i)
        el <- as.integer(gsub("[[:alpha:]]", "", i))
        nd[, la] ^ el
      })
      po <- as.data.frame(do.call("cbind", l))
      names(po) <- elev
      nd <- cbind(nd, po)
    }
    # add interactions
    if (length(vi)) {
      m <- lapply(vi, function(i) {
        va <- strsplit(i, "\\.")[[1]]
        nd[, va[1]] * nd[, va[2]]
      })
      int <- as.data.frame(do.call("cbind", m))
      names(int) <- vi
      nd <- cbind(nd, int)
    }
    nd <- nd[, names(x), drop = FALSE] # reorder according to model.frame
    if (length(vars) == 2) {
      res <- split(nd, nd[, vars[2]]) # returns list of lists if interaction
      names(res) <- 1:n2
      res
    } else{
      list(nd) # a length 1 list if simple
    }
  })
  names(ans) <- sapply(varset, paste, collapse = ".")
  ans
}

###################################################################
### Range hull function
# Studies relation between nb of locations and hull size of occurrences

rangeHull<-function(x,species="species",breaks=50,plot=TRUE){
  n<-1:nrow(x)
  n<-round(seq(min(n),max(n),length.out=breaks),0)
  l<-lapply(n,function(i){
    print(i)
    xx<-x[sample(1:nrow(x),i,replace=FALSE),]
    st_union(xx)
  })

#obs<-rev(sort(table(x$recordedBy)))[1:10]

  m<-do.call("c",l)
  m<-st_convex_hull(m)
  
  oldpar<-par()
  on.exit(suppressWarnings(par(oldpar)))
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
    list(hull=m,nobs=nrow(x),a=coef(model)["a"],reach=a/coef(model)["a"],n=n,model=model)
  }else{
    list(hull=st_as_sf(m),nobs=nrow(x),a=NA,reach=NA,n=n,model=model)
  }
}

### plots asymptotic hull
#rh<-rangeHull(occs,breaks=100)

stabHull<-function(species,path="/data/sdm_rbq/stab",rh=NULL){
  areas<-st_area(rh$hull)
  expo<-floor(log(max(as.numeric(areas)),base=10))
  png(file.path(path,paste0(gsub(" ","_",species),"_stab.png")),width=9,height=8,units="in",res=200,pointsize=11)
  par(mar=c(4,4.5,1,1))
  ylim<-c(0,max(c(rh$a,max(areas)))*1.05)
  plot(rh$n,areas,ylim=ylim,yaxt="n",xaxt="n",xlab="",ylab="",col="black",bg=adjustcolor("black",0.25),pch=21,cex=2.5)
  axis(1,tcl=-0.3,mgp=c(2,0.5,0),cex.axis=1.5,lwd=0,lwd.ticks=1)
  axis(2,las=2,tcl=-0.3,mgp=c(2,0.5,0),cex.axis=1.5,lwd=0,lwd.ticks=1,at=pretty(areas),labels=pretty(areas)/10^expo)
  abline(h=rh$a,lty=3,lwd=8,col="black")
  vs<-0:max(rh$n)
  lines(vs,predict(rh$model,data.frame(X=vs)),lwd=8,col="black")
  box(col="snow2",lwd=5)
  mtext(side=1,line=2.5,cex=2.5,font=2,text="Number of locations sampled")
  mtext(side=2,line=2,cex=2.5,font=2,text=bquote(bold("Convex Hull Area "%*%~10^.(expo)~km^2)))
  mtext(side=1,line=-5,adj=0.8,cex=4.5,font=2,text=round(rh$reach,2),col="grey20")
  #text(par("usr")[1],ylim[2]*1.02,label="Asymptote",adj=c(-0.1,1),cex=2.5,font=2)
  dev.off()
}

#stabHull(species,rh=rh)

region_mesh<-function(mesh){
  x<-unique(as.vector(mesh$segm$int$idx))
  reg<-st_polygon(list(Mesh$loc[c(x,1),]))
  st_as_sf(st_sfc(reg),crs=st_crs(na))
}



#plot(Mesh,asp=TRUE)
#lines(Mesh$loc[x,1],Mesh$loc[x,2],lwd=5,col="red")
#meshregion<-st_polygon(list(Mesh$loc[c(x,1),]))
#meshregion<-st_as_sf(st_sfc(meshregion),crs=st_crs(na))
#plot(st_geometry(meshregion),add=TRUE,border="red",lwd=5)

#############################################################
### Mesh to sp function #####################################

inla.mesh2sp <- function(mesh) {
#crs <- inla.CRS(inla.CRSargs(mesh$crs))
crs <- inla.CRS(mesh$crs)
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


### objects

library(FRutils)

colmean<-c("#CCCCCC","#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00","darkorange")
colmean<-colo.scale(1:200,colmean)

naplot<-function(lwd=0.25){
  plot(st_geometry(na),border="white",lwd=lwd,axes=FALSE,add=TRUE)
  plot(st_geometry(nalakes),col="white",border=adjustcolor("black",0.15),lwd=lwd,axes=FALSE,add=TRUE)
  plot(st_geometry(coast),lwd=lwd,border=gray(0,0.45),axes=FALSE,add=TRUE)
}

getobs<-function(sp){
  be<-unlist(ed[match(sp,ed$scientific_name),c("start","end")])
  ### all obs for period
  obs<-d[md>=be[1] & md<=be[2],]
  ### remove duplicate obs
  obs<-unique(obs,by=c("recordedBy","species","cell"))
  ### species obs
  spobs<-obs[species==sp,]
  nndist<-knn.dist(as.matrix(spobs[,c("x","y")]),k=1)[,1]
  spobs<-spobs[nndist<=1200,]
  st_as_sf(spobs,coords=c("x","y"),crs=crsr)
}


#xx<-getobs("Setophaga petechia")

overlapHist<-function(x="I",sp=NULL,th=0.85,n=0){
  df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
  if(is.null(sp)){
    path<-"/data/sdm_rbq/graphics/overlap.png"
    m<-NA
    bg<-"white"
  }else{
    path<-file.path("/data/sdm_rbq/overlap",paste0(gsub(" ","_",sp),"_over.png"))
    bg<-"white"#"lemonchiffon2"
  }
  png(path,width=12,height=9,units="in",res=400)
  par(bg=bg)
  ag<-df[,x]
  ag<-ag[which(df$reach>=th & df$n>=n)]
  brks<-seq(0.0,1,by=0.1)
  h<-hist(ag,xlim=c(0,1),breaks=brks,border="white",lwd=5,col=adjustcolor("seagreen",0.5),main="",xaxt="n",yaxt="n",xlab="",ylab="")
  yat<-pretty(h$counts)
  yat<-c(yat,yat[length(yat)]+diff(yat)[1])
  yat<-unique(c(0,yat))
  lab<-yat
  lab[1]<-""
  axis(1,at=brks,mgp=c(3,1,0),pos=c(0,0),col="grey50",col.axis="black",col.ticks="grey50",lwd=5,cex.axis=2,font.axis=2,tcl=0.5)
  if(is.null(sp)){
    axis(2,at=yat,labels=lab,mgp=c(2.25,0.75,0),pos=c(-0.0,0),col="grey50",col.axis="black",col.ticks="grey50",lwd=5,cex.axis=2,font.axis=2,las=2,tcl=0.5,xpd=FALSE)
    axis(3,at=brks,labels=rep("",length(brks)),mgp=c(2,1,0),pos=c(par("usr")[4],0),col="grey50",col.axis="black",col.ticks="grey50",lwd=5,cex.axis=1.5,font.axis=2,las=2,tcl=0.5)
    axis(4,at=yat,labels=rep("",length(yat)),mgp=c(2,1,0),pos=c(1,0),col="grey50",col.axis="black",col.ticks="grey50",lwd=5,cex.axis=1.5,font.axis=2,las=2,tcl=0.5)
    mtext(side=2,line=1,text="Number of species",font=2,cex=3)
  }
  if(!is.null(sp)){
    m<-match(sp,df$species)
    lines(rep(df[,x][m],2),c(0,par("usr")[4]),lwd=8,col="black")
  }else{
    lines(rep(median(ag,na.rm=TRUE),2),c(0,par("usr")[4]),lwd=8,col="black")
    legend("topleft",legend="Median",,lwd=8,col="black",bty="n",inset=c(0.1,0.1),cex=3,text.font=2)
  }
  mtext(side=1,line=2.5,text="Overlap (I)",font=2,cex=3)
  
  dev.off()
  list(m=m,path=path)
}


getCC0links<-function(species,license=c("cc0")){  
  sp<-gsub(" ","%20",species)
  cc<-paste(license,collapse="0%2C")
  urlsearch<-paste0("https://api.inaturalist.org/v1/taxa?q=",sp,"&order=desc&order_by=observations_count")
  x<-fromJSON(urlsearch)$results
  m<-match(species,x$name)
  if(is.na(m)){
    m<-match(species,x$matched_term)
    if(is.na(m)){
      m<-grep(species,x$matched_term)[1]
    }  
  }
  if(is.na(m)){
    warning(sprintf("No match for %s",species))
    taxonid<-x$id[1]
  }else{
    taxonid<-x$id[m]
  }
  x<-fromJSON(paste0("https://api.inaturalist.org/v1/observations?photo_license=",cc,"&taxon_id=",taxonid,"&&quality_grade=research&order=desc&order_by=created_at"))
  if(x$total_results==0){
    return(NA)
  }else{
    x<-x$results
  }
  users<-x$user[,c("login","name")]
  pics<-do.call("rbind",lapply(seq_along(x$observation_photos),function(i){
    res1<-x$observation_photos[[i]]$photo[,c("url","license_code","attribution")]
    res2<-x$observation_photos[[i]]$photo$original_dimensions[,c("width","height")]
    res<-cbind(res1,res2)
    #res<-res[which(res$width>205 & res$height>205),]
    #if(nrow(res)>0){
    res<-res[1,] # keep first one
    #}
    cbind(id=x$id[i],res,users[rep(i,nrow(res)),])
  })) 
  pics$url<-gsub("/square","/medium",pics$url)
  pics<-pics[which(pics$width>205 & pics$height>205),]
  pics
}
  
#getCC0links("Melinis repens") 

roundIm<-function(url,file,path,link=NULL,license=NULL,author=NULL,open=FALSE){
  if(url%in%c("",NA)){
    im<-image_blank(width=500, height=500, color = "")
  }else{
    im <- image_read(url)
  }
  wh<-c(image_info(im)$width,image_info(im)$height)
  mi<- min(wh)
  wm<-which.max(wh)
  if(wh[1]!=wh[2]){
    if(wm==1){
      geom<-paste0(mi, "x", mi, "+",abs(diff(wh))/2,"+","0")
    }else{
      geom<-paste0(mi, "x", mi, "+","0","+",abs(diff(wh))/2)
    }
    im<-image_crop(im, geometry=geom,repage=TRUE)
  }
  wh<-c(image_info(im)$width,image_info(im)$height)
  mask <- image_draw(image_blank(wh[1],wh[2]))
  symbols(wh[1]/2,wh[2]/2,circles=(wh[1]/2)-2,bg="#000000",inches=FALSE,add=TRUE)
  dev.off()
  organism<-image_composite(im,mask, operator="copyopacity") 
  organim<-image_trim(organism)
  h<-image_info(organism)$height
  organism<-image_border(organism,"grey95","50x50",operator="over")
  #organism<-image_draw(organism)
  #symbols(image_info(organism)$height/2,image_info(organism)$width/2,circles=h/2,bg="transparent",fg="grey99",inches=FALSE,lwd=20,add=TRUE)
  #dev.off()
  #organism<-image_trim(organism,fuzz=0)
  #organism<-image_border(organism,fuzz=0)
  organism<-image_fill(organism,"none",point="+5+5",fuzz=1)
  organism<-image_trim(organism,fuzz=0)
  organism<-image_scale(organism,"250")
  image_write(organism,file.path(path,paste0(file,"_pic.png")),format="png",comment=link)
  tags<-c("-overwrite_original",paste0("-ImageDescription=",link),paste0("-Copyright=",license))
  #print(tags)
  exif_call(path=file.path(path,paste0(file,"_pic.png")),args=tags)
}

#lapply(seq_along(pics),function(i){
#  roundIm(url=pics[[i]],file=names(pics)[i],path=path,open=FALSE)
#})

#spp<-
#speciesPics<-function(species=spp)

coloScale<-function(
    vals,#values(r)[,1],
    xy=c(1.015,1),
    cols=rev(terrain.colors(200)),
    at=NULL,
    labels=NULL,
    width=0.015,
    height=1,
    ticks=c(0.15,0),
    labs=c(0.25),
    cex=0.9,
    lim=FALSE
  ){
  vals<-range(vals,na.rm=TRUE)
  usr<-par("usr")  
  xleft<-usr[1]+diff(usr[1:2])*xy[1]
  ytop<-usr[3]+diff(usr[3:4])*xy[2] 
  xright<-xleft+(diff(usr[1:2])*width)
  ybottom<-ytop-(diff(usr[3:4])*height)
  rect(xleft=xleft,ybottom=ybottom,xright=xright,ytop=ytop,xpd=TRUE,border=NA,lwd=1,col=NA)
  yy<-seq(ybottom,ytop,length.out=length(cols)+1)
  for(i in seq_along(cols)){
    rect(xleft=xleft,ybottom=yy[i+1],xright=xright,ytop=yy[i],xpd=TRUE,border=NA,col=cols[i])  
  }
  
  if(is.null(at)){
    at<-pretty(vals)
    at<-at[at>=min(vals) & at<=max(vals)]
    if(lim){
      at<-sort(unique(c(vals,at)))  
    }
    labels<-at
  }
  
  at2<-scales::rescale(c(at,vals),to=c(ybottom,ytop))[1:length(at)]
  
  offset<-(xright-xleft)*ticks*c(-1,1)
  offset2<-(xright-xleft)*labs
  for(i in seq_along(at2)){
    #lines(rep(xleft,2)+offset,rep(at2[i],2),col="black",xpd=TRUE)  
    lines(rep(xright,2)+offset,rep(at2[i],2),col="black",xpd=TRUE) 
    text(xright+offset2,at2[i],label=labels[i],cex=cex,adj=c(0,0.5),xpd=TRUE) 
  }
  
  
}


map<-function(name, obs=TRUE, cols=TRUE, trans=NULL, ...){
  opts<-par("mfrow")
  par(mfrow=n2mfrow(length(name),asp=2))
  on.exit(par(opts))
  #name<-deparse(substitute(name))
  for(i in name){
  plg<-list(cex=1.5,shrink=0.9)
  mar<-c(0,0,0,0)
  col<-if(cols){colmean}else{rev(terrain.colors(200))}
  if(i%in%names(sdm)){
    if(i=="spacemean"){
        vals<-unlist(global(sdm[[i]],"range",na.rm=TRUE))
        col<-colo.scale(seq(vals[1],vals[2],length.out=200),c("grey20","blue4","dodgerblue",colmean[1],"tomato2","red4","grey20"),center=TRUE)
    }
    plot(mask(mask(sdm[[i]],region),na),col=col,plg=plg,mar=mar,...)
  }else{
    if(i=="pred"){
      plot(mask(mask(exp(sdm[["linkmean"]]-sdm[["spacemean"]]),region),na),col=col,plg=plg,mar=mar,...) 
    }else{
      if(i=="ebird"){
        ebirdpath<-file.path("/data/predictors_sdm/expert_maps/eBird/abundance",paste0(gsub(" ","_",sp),"_ebird2.tif"))
        ebird<-rast(ebirdpath)
        if(!is.null(trans)){
          ebird<-trans(ebird)
        }
        plot(mask(mask(ebird,vect(region)),na),col=col,plg=plg,mar=mar,...) 
      }else{
        plot(mask(mask(r[[i]],region),na),col=col,plg=plg,mar=mar,...) 
      }
    }
  }
  naplot(lwd=1)
  if(obs){
    plot(st_geometry(occs),add=TRUE)
  }
}
}

#map(c("mean","ebird"),obs=FALSE)





##########################################################
### plot marginale effects ###############################

plot_marginal_effects<-function(m,sp,dmesh,nsims=100,open=FALSE){

#m<-mpp
class(m)<-"inla"
samples<-inla.posterior.sample(nsims,m,num.threads="2:2")

params<-dimnames(m$model.matrix)[[2]]
nparams<-sapply(params,function(i){
  match(paste0(i,":1"),row.names(samples[[1]]$latent)) 
}) 
tab<-table(sapply(strsplit(row.names(samples[[1]]$latent),":"),"[",1))
tab
nweights<-grep("i:",row.names(samples[[1]]$latent))

fun<-function(i){ # gets median from where the species is found
  median(i[dmesh$spobs>0]) 
  #mean(i)
}

mm<-as.data.frame(as.matrix(m$model.matrix)[,-1,drop=FALSE])[1:(tab[1]/2),,drop=FALSE] # what's in the model.matrix is related to Predictor, not APredictor
mm<-as.data.frame(dmeshPred[,names(mm)]) # I think better to use this
nd<-newdata2(mm,n=100,n2=5,fun=mean)
#print(dmesh$spobs)
v<-names(nd)
nd<-unlist(nd,recursive=FALSE)

#mm<-as.data.frame(dmeshPred[,c("tmean","grass_esa")]) # I think better to use this
#nd<-newdata2(mm,n=100,n2=5,fun=fun)
#v<-names(nd)
#nd<-unlist(nd,recursive=FALSE)

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
  p1<-base::t(apply(p1,1,function(i){c(quantile(i,0.025,na.rm=TRUE),mean(i,na.rm=TRUE),quantile(i,0.975,na.rm=TRUE))}))
  p2<-base::t(apply(p2,1,function(i){c(quantile(i,0.025,na.rm=TRUE),mean(i,na.rm=TRUE),quantile(i,0.975,na.rm=TRUE))}))
  p1<-exp(p1)
  p2<-exp(p2)
  p1<-setNames(as.data.frame(p1),c("low","mean","high"))
  add<-strsplit(names(nd)[k],"\\.")[[1]]
  data.frame(vars=add[1],block=add[2],p1)
})
names(marginals)<-names(nd)


preds<-cbind(do.call("rbind",marginals),do.call("rbind",nd))


ma<-match(names(preds),names(means))
ma<-ma[!is.na(ma)]
preds[names(means)[ma]]<-lapply(names(means)[ma],function(i){
  backScale(preds[[i]],i)
})


preds<-split(preds,preds$vars)
preds<-lapply(preds,function(i){
  if(all(is.na(i$block))){
    list(i)
  }else{
    split(i,i$block)
  }
})


pathg<-paste0("/data/sdm_rbq/marginaleffects/",gsub(" ","_",sp),"_me.png")
mfrow<-n2mfrow(length(preds)+2,asp=1.5)
png(pathg,width=3.5*mfrow[2],height=3.5*mfrow[1],units="in",res=300)
par(mfrow=mfrow,mar=c(2.5,1.25,0.5,0.5),oma=c(1,2,5,1))
lapply(names(preds),function(i){
  xvar<-strsplit(i,"\\.")[[1]][1]
  xvar2<-strsplit(i,"\\.")[[1]][2]
  xlim<-range(preds[[i]][[1]][,xvar])
  ywide<-unlist(preds)
  ylim<-range(as.numeric(ywide[grep("\\.mean",names(ywide))]))*c(0.5,1.5)
  #ylim<-range(unlist(do.call("rbind",preds[[i]])[,c("low","mean","high"),drop=FALSE]))
  xaxt<-if(i=="logdistance"){"n"}else{"s"}
  yaxt<-ifelse(match(i,names(preds))%in%seq(1,mfrow[1]*mfrow[2],by=mfrow[2]),"s","n")
  plot(0.1,0.1,xlab="",ylab="",xlim=xlim,ylim=ylim,type="n",mgp=c(2,0.15,0),tcl=-0.1,log="",xaxt=xaxt,yaxt=yaxt)
  box(col="white",lwd=2)
  grid(lty=3,col=adjustcolor("black",0.2))
  if(i=="logdistance"){
    axis(1,at=pretty(xlim,20),labels=round(10^pretty(xlim,20)-1,1),mgp=c(2,0.15,0),tcl=-0.1,cex.axis=0.5)
  }
  ### density of values where species is observed
  vals<-backScale(dmeshPred[,i],i)
  brks<-seq(min(vals),max(vals),length.out=20)

  h<-hist(vals,breaks=brks,plot=FALSE)
  h$density<-(h$density/max(h$density))
  h$density<-h$density/(max(h$density)/ylim[2])
  lapply(seq_along(h$density),function(j){
    rect(xleft=h$breaks[j],ybottom=0,xright=h$breaks[j+1],ytop=h$density[j],border=NA,col=adjustcolor("seagreen",0.25))
  })
  
  h<-hist(vals[dmesh$spobs>0],breaks=brks,plot=FALSE)
  h$density<-(h$density/max(h$density))
  h$density<-h$density/(max(h$density)/ylim[2])
  lapply(seq_along(h$density),function(j){
    rect(xleft=h$breaks[j],ybottom=0,xright=h$breaks[j+1],ytop=h$density[j],border=NA,col=adjustcolor("orange",0.30))
  })

  lapply(seq_along(preds[[i]]),function(jj){
    j<-preds[[i]][[jj]]
    x<-j[,xvar]
    polygon(c(x,rev(x),x[1]),c(j[,"low"],rev(j[,"high"]),j[,"low"][1]),col=adjustcolor("black",0.13),border=NA)
    lines(x,j[,"mean"],lwd=3,col=adjustcolor("black",0.5))
    lines(x,j[,"low"],lwd=1,lty=1,col=adjustcolor("black",0.25))
    lines(x,j[,"high"],lwd=1,lty=1,col=adjustcolor("black",0.25))
    mtext(side=1,line=1.5,text=xvar,font=2,cex=1.25)
    #box(col="grey70")
  })
  if(length(preds[[i]])>1){
    temp<-do.call("rbind",preds[[i]])
    var2<-strsplit(temp$vars[1],"\\.")[[1]][2]
    legend("top",title=var2,legend=round(as.numeric(unique(temp[,var2])),2),lwd=1,col=seq_along(preds[[i]]),bty="n",cex=1.25)
  }
})
posterior<-function(param,name){
  i<-grep(param,names(m$marginals.hyperpar),ignore.case=TRUE)
  xy<-m$marginals.hyperpar[[i]]
  xy2<-xy[rev(1:nrow(xy)),]
  xy2[,2]<-0
  xy2<-rbind(xy,xy2)
  plot(xy[,1],xy[,2],xlim=c(min(xy[,1]),xy[,1][max(which(xy[,2]>(max(xy[,2])*0.025)))+1]),type="n",xlab="",ylab="",mgp=c(2,0.25,0),tcl=-0.1)
  box(col="white",lwd=2)
  grid(lty=3,col=adjustcolor("black",0.2))
  polygon(xy2[,1],xy2[,2],border=NA,col=adjustcolor("black",0.10))
  lines(xy[,1],xy[,2],lwd=1,col=adjustcolor("black",0.25))
  abline(v=m$summary.hyperpar[i,][,"mean"],lwd=3,col=adjustcolor("black",0.50))
  mtext(side=1,line=2,text=name,font=2,cex=1.25)
  #box(col="grey70")
}
posterior("Range","Range (km)")
posterior("stdev","SD")

par(new=TRUE,mfrow=c(1,1),mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(0,0,axes=FALSE,type="n",bty="n")
legend("top",ncol=2,legend=c("Mean (predicted or parameter)","CI (95%) or posterior distribution","Density of values in study area","Density of values where species observed"),pch=c(NA,22,22,22),pt.cex=c(NA,2,2,2),lwd=c(3,NA,NA,NA),pt.lwd=c(3,NA,NA,1),inset=c(0,0),col=c(adjustcolor("black",0.50),adjustcolor("black",0.25),adjustcolor("seagreen",0.25),adjustcolor("orange",0.25)),pt.bg=c(NA,adjustcolor("black",0.15),adjustcolor("seagreen",0.25),adjustcolor("orange",0.25)),cex=1.25,bty="n")

dev.off()
if(open){
  system(paste("code",pathg), wait=FALSE)
}
}

#plot_marginal_effects(mpp,open=TRUE)

#########################################################
### pixelation ##########################################

pixelate<-function(sdm,sp){
  png(paste0("/data/sdm_rbq/figures/",gsub(" ","_",sp),"_pixelation.png"),width=5,height=5,res=400,units="in")
  sdm2<-sdm[[grep("sample",names(sdm))]]
  #sdm2<-aggregate(sdm2,2,fun=function(i){sample(i,1)})
  a<-as.array(sdm2)
  pixels<-sapply(unlist(asplit(asplit(a,1:2),1),recursive=FALSE,use.names=FALSE),function(i){
    #i<-i[i>=quantile(i,0.025,na.rm=TRUE) & i<=quantile(i,0.975,na.rm=TRUE)]
    sample(i,size=1)
  })
  pixelation<-setValues(sdm2[[1]],pixels)
  pixelation<-resample(pixelation,sdm[[1]])
  xy=c(0.91,0.42)
  height=0.37
  cex=0.75
  plot(pixelation,col=colmean,mar=c(0,0,0,0),legend=FALSE,axes=FALSE)
  coloScale(vals=range(values(pixelation),na.rm=TRUE),cols=colmean,xy=xy,height=height,cex=cex)
  naplot()
  dev.off()
}

#pixelate(sdm)
#########################################################
### animation ###########################################

animate<-function(sdm,sp){
  samps<-sample(grep("sample",names(sdm),value=TRUE),40)
  zlim<-range(values(sdm[[samps]]),na.rm=TRUE)
  xy=c(0.91,0.42)
  height=0.37
  cex=0.75
  ims<-lapply(seq_along(samps),function(i){
    png(paste0("/data/sdm_rbq/temp/animation",sp,i,".png"),width=5,height=5,res=400,units="in")
    #im <- image_graph(width = 900, height = 800, res = 96)
    #map(samps[i],range=zlim,obs=FALSE)
    plot(sdm[[samps[i]]],col=colmean,mar=c(0,0,0,0),legend=FALSE,axes=FALSE)
    naplot()
    coloScale(vals=zlim,cols=colmean,xy=xy,height=height,cex=cex)
    dev.off()
  })
  lf<-list.files("/data/sdm_rbq/temp",pattern=paste0("animation",sp),full=TRUE)
  img<-lapply(lf,image_read)
  img<-do.call("c",img)
  img<-image_scale(img,"x450")
  animation <- image_animate(img, fps = 10, optimize = TRUE)
  image_write(animation,paste0("/data/sdm_rbq/figures/",gsub(" ","_",sp),"_gif.gif"))
}

#animate(sdm)

