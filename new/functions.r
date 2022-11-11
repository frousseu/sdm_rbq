
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
    ), "_"))) # all unique variable names
  vs <-
    unique(unlist(strsplit(names(x), "_"))) # all vars or components of interactions
  vi <- names(x)[grep("_", names(x))] # all interactions
  fix <-
    as.data.frame(as.list(apply(x[, vn, drop = FALSE], 2, FUN = fun))) # mean values of each var
  # build variables set (corresponds to all graphs)
  temp <- strsplit(names(x), "_")
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
        va <- strsplit(i, "_")[[1]]
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
  names(ans) <- sapply(varset, paste, collapse = "_")
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