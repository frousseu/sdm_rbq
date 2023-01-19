
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

naplot<-function(){
  plot(st_geometry(na),border="white",lwd=0.25,axes=FALSE,add=TRUE)
  plot(st_geometry(nalakes),col="white",border=adjustcolor("black",0.15),lwd=0.25,axes=FALSE,add=TRUE)
  plot(st_geometry(coast),lwd=0.25,border=gray(0,0.45),axes=FALSE,add=TRUE)
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
    vals=values(r)[,1],
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
