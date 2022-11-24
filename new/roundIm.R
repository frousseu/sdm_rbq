
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
  plot(0:1,0:1,type="n",xaxt="n",yaxt="n",bty="n",bg="white")
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



