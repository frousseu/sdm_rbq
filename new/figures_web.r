
# this is the part for creating figures that comes from the old visualizations.r file

# find -mtime -1 -ls | wc -l

# nohup Rscript figures_web.r --no-save > verbose.out 2>&1 &


library(tidyr)
library(foreach)
library(doParallel)
library(FRutils)
library(berryFunctions)
library(magick)
library(jsonlite)
library(exiftoolr)
library(colorspace)

options(vsc.dev.args = list(width = 1000, height = 800))

code<-readLines("/data/sdm_rbq/models.r",n=360) # load everything first
code<-code[1:grep("Data loaded",code)]
code<-paste(code,collapse="\n")
source(textConnection(code))


########################################################################
### Plot global and individual overlap #################################
########################################################################

df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")

overlapHist(x="I",th=0.85,n=50)
#o<-overlapHist(sp="Columbia livia")
#plot(image_read(o$path),"#FFFFFF00")

cl<-makeCluster(20)
registerDoParallel(cl)
foreach(i=df$species,.packages=c("sf")) %dopar% {
  overlapHist(sp=i)
}
stopCluster(cl)



######################################################################
### Species pics #####################################################
######################################################################

if(FALSE){

path<-"/data/sdm_rbq/pics"

df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
#sps<-c("Buteo jamaicensis","Buteo lagopus","Buteo platypterus")
sps<-sort(df$species)#[1:18]
#sps<-sps[544:length(sps)]
#i<-c("Acanthis hornemanni")
#sps<-i#c("Columba livia")

 # run if needed
lapply(sps,function(i){
  print(i)
  links<-getCC0links(i,license=c("cc0"))
  #Sys.sleep(0.05)
  if(identical(NA,links)){
    url<-"https://inaturalist-open-data.s3.amazonaws.com/photos/246765575/large.jpg"
    url<-""
    link<-"https://www.inaturalist.org/observations/143844420"
    link<-""
    author<-""
    license<-""

  }else{
    url<-links$url[1]
    link<-paste0("https://www.inaturalist.org/observations/",links$id[1])
    if(links$license_code[1]=="cc0"){
      if(links$name[1]%in%c("",NA)){
        author<-links$login[1]
      }else{
        author<-links$name[1]
      }
      license<-paste0("(c) ",author,", ",links$attribution[1]," (CC0)")
    }else{
      license<-links$attribution[1]
    } 
    
  }
  roundIm(url=url,file=gsub(" ","_",i),path=path,link=link,license=license,open=FALSE)
  print(i)
})

}

#exif_read(path = "/data/sdm_rbq/pics/Acanthis_hornemanni_pic.png", #tags = "*Copyright*")$Copyright


###################################################################
###################################################################
### Stabilization graphs ##########################################
###################################################################
###################################################################

if(FALSE){ # already ran with modelss

lf<-list.files("/data/sdm_rbq/stab/")
lf<-sapply(strsplit(lf,"_"),function(i){paste(i[1:2],collapse=" ")})
df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
#df<-df[df$n>=5,]
#sp<-sort(df$species[!df$species%in%lf])
sp<-df$species
#i<-"Buteo jamaicensis"

#cl<-makeCluster(2)
#registerDoParallel(cl)
#foreach(i=sp) %dopar% {
 
lapply(sp,function(i){
  print(i)
  sPoints<-getobs(i)
  rh<-rangeHull(sPoints,species=i,breaks=200)
  stabHull(i,rh=rh)
})

}
#}
#stopCluster(cl)


##########################################################################
### Web maps and side by side ############################################
##########################################################################

#list.files(tools::R_user_dir("ebirdst"),recursive=TRUE)
#options(vsc.dev.args = list(width = 1500, height = 800))
#colmean<-c("#CCCCCC","#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00","darkorange")
#colmean<-colo.scale(1:200,colmean)
#coast<-st_union(na)

#naplot<-function(){
#  plot(st_geometry(na),border="white",lwd=0.25,add=TRUE)
#  plot(st_geometry(nalakes),col="white",border=adjustcolor("black",0.15),lwd=0.25,add=TRUE)
#  plot(st_geometry(coast),lwd=0.25,border=gray(0,0.45),add=TRUE)
#}


df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
spc<-df$species[which(df$I>=0.75 & df$I<=0.95 & df$reach>=0.85 & df$n>=100)]
spc
#df<-df[df$date>="2023-02-12 00:00:00",]

lsp<-unique(df$species)#[1]
#i<-"Catharus bicknelli"
#lsp<-i


#plan(multicore,workers=5)
#cores<-nbrOfWorkers() # get nbr of workers from the chosen plan
#chunks <- split(1:length(lsp), rep(1:cores, each=ceiling(length(lsp)/cores))[1:length(lsp)])

cl<-makeCluster(20)
registerDoParallel(cl)
foreach(i=lsp,.packages=c("sf","terra","magick","FRutils","data.table","FNN","colorspace")) %dopar% {
#res<-future_lapply(chunks,function(chunksi){
  print(i)
  sp<-i
  sdm<-rast(file.path("/data/sdm_rbq/rasters",paste0(gsub(" ","_",sp),"_birds.tif")))[["mean"]]
  linkmean<-rast(file.path("/data/sdm_rbq/rasters",paste0(gsub(" ","_",sp),"_birds.tif")))[["linkmean"]]
  linksd<-rast(file.path("/data/sdm_rbq/rasters",paste0(gsub(" ","_",sp),"_birds.tif")))[["linksd"]]
  spacemean<-rast(file.path("/data/sdm_rbq/rasters",paste0(gsub(" ","_",sp),"_birds.tif")))[["spacemean"]]
  predmean<-exp(linkmean-spacemean)
  ebirdpath<-file.path("/data/predictors_sdm/expert_maps/eBird/abundance",paste0(gsub(" ","_",sp),"_ebird2.tif"))
  if(file.exists(ebirdpath)){  
    ebird<-rast(ebirdpath)
  }else{
    return(NULL) #next
  }

  fname<-gsub(" ","_",sp)
  colspacemean<-colo.scale(seq(min(values(spacemean),na.rm=TRUE),max(values(spacemean),na.rm=TRUE),length.out=200),c("grey20","blue4","dodgerblue",colmean[1],"tomato2","red4","grey20"),center=TRUE)
   colssd<-rev(sequential_hcl(5, palette = "Greens"))
   colssd[1]<-colmean[1]
   colssd<-colo.scale(1:200,c(colssd,"black"))

  xy=c(0.91,0.42)
  height=0.37
  cex=0.75


  png(file.path("/data/sdm_rbq/temp",paste0(fname,"_obs.png")),width=5,height=5,res=400,units="in")
  plot((sdm)^(1/1),col=colmean[1],mar=c(0,0,0,0),legend=FALSE,axes=FALSE)
  naplot()
  plot(st_geometry(getobs(sp)),pch=16,cex=0.9,col=adjustcolor(colmean[95],0.5),add=TRUE)
  dev.off()
  png(file.path("/data/sdm_rbq/temp",paste0(fname,"_sdm.png")),width=5,height=5,res=400,units="in")
  plot((sdm)^(1/1),col=colmean,mar=c(0,0,0,0),legend=FALSE,axes=FALSE)
  naplot()
  vals<-values(sdm)
  coloScale(vals=vals,cols=colmean,xy=xy,height=height,cex=cex)
  dev.off()
  expo<-1/2
  ebird2<-aggregate(ebird,9,na.rm=TRUE)^expo
  #ebird2<-ebird;expo<-1/1
  png(file.path("/data/sdm_rbq/temp",paste0(fname,"_ebird.png")),width=5,height=5,res=400,units="in")
  plot(ebird2,col=colmean,mar=c(0,0,0,0),legend=FALSE,axes=FALSE)
  naplot()
  vals<-range(values(ebird2),na.rm=TRUE)^expo
  at<-pretty(vals)
  at<-at[at>=min(vals) & at<=max(vals)]
  #at<-unique(sort(c(vals,at)))
  labels<-round(at^(1/expo),1)
  coloScale(vals=vals,cols=colmean,xy=xy,height=height,cex=cex,at=at,labels=labels)
  dev.off()
  png(file.path("/data/sdm_rbq/temp",paste0(fname,"_sd.png")),width=5,height=5,res=400,units="in")
  plot((linksd)^(1/1),col=colssd,mar=c(0,0,0,0),legend=FALSE,axes=FALSE)
  naplot()
  vals<-values(linksd)
  coloScale(vals=vals,cols=colssd,xy=xy,height=height,cex=cex)
  dev.off()
  png(file.path("/data/sdm_rbq/temp",paste0(fname,"_pred.png")),width=5,height=5,res=400,units="in")
  plot((predmean)^(1/1),col=colmean,mar=c(0,0,0,0),legend=FALSE,axes=FALSE)
  naplot()
  vals<-values(predmean)
  coloScale(vals=vals,cols=colmean,xy=xy,height=height,cex=cex)
  dev.off()
  png(file.path("/data/sdm_rbq/temp",paste0(fname,"_spatial.png")),width=5,height=5,res=400,units="in")
  plot((spacemean)^(1/1),col=colspacemean,mar=c(0,0,0,0),legend=FALSE,axes=FALSE)
  naplot()
  vals<-values(spacemean)
  coloScale(vals=vals,cols=colspacemean,xy=xy,height=height,cex=cex)
  dev.off()


  if(FALSE){
  im<-image_trim(image_read(file.path("/data/sdm_rbq/temp",paste0(fname,"_obs.png"))))
  im<-image_border(im,"white","1600")
  im2<-image_trim(image_read(file.path("/data/sdm_rbq/temp",paste0(fname,"_sdm.png"))))
  im2<-image_fill(im2,"#FFFFFF00","+5+5")
  im<-image_composite(im,im2,gravity="east",offset="+0+0")
  im<-image_trim(im)
  im<-image_border(im,"white","1600")
  im3<-image_trim(image_read(file.path("/data/sdm_rbq/temp",paste0(fname,"_ebird.png"))))
  im3<-image_fill(im3,"#FFFFFF00","+5+5")
  im<-image_composite(im,im3,gravity="east",offset="+0+0")
  im<-image_trim(im)
  im<-image_border(im,"white","1600")
  im4<-image_trim(image_read(file.path("/data/sdm_rbq/temp",paste0(fname,"_sd.png"))))
  im4<-image_fill(im4,"#FFFFFF00","+5+5")
  im<-image_composite(im,im4,gravity="east",offset="+0+0")
  im<-image_trim(im)
  im<-image_border(im,"white","1600")
  im5<-image_trim(image_read(file.path("/data/sdm_rbq/temp",paste0(fname,"_spatial.png"))))
  im5<-image_fill(im5,"#FFFFFF00","+5+5")
  im<-image_composite(im,im5,gravity="east",offset="+0+0")
  im<-image_trim(im)

  im<-image_border(im,"white","50x50")
  im<-image_annotate(im,sp,size=120,color="black",gravity="northwest",location="+10+10",boxcolor="#FFFFFFAA",degrees=0,weight=700)
  m<-match(sp,df$species)
  im<-image_annotate(im,paste0("n = ",df$n[m]),size=80,color="black",gravity="west",location="+10-100",degrees=0,weight=700)
  im<-image_annotate(im,paste0("reach = ",round(df$reach,2)[m]),size=80,color="black",gravity="west",location="+10+0",degrees=0,weight=700)
  im<-image_annotate(im,paste0("cor = ",round(df$pearson[m],2)),size=80,color="black",gravity="west",location="+10+100",degrees=0,weight=700)
  im<-image_annotate(im,paste0("I = ",round(df$I[m],2)),size=80,color="black",gravity="west",location="+10+200",degrees=0,weight=700)
  im<-image_annotate(im,"mapSpecies",size=100,color="#00000044",gravity="northwest",location="+750+300",boxcolor="#FFFFFFAA",degrees=0,weight=700)
  im<-image_annotate(im,"eBird",size=100,color="#00000044",gravity="northwest",location="+2550+300",boxcolor="#FFFFFFAA",degrees=0,weight=700)
  im<-image_scale(im,"x500")
  image_write(im,file.path("/data/sdm_rbq/comparison",paste0(fname,"_comparison.png")))
  }
}
stopCluster(cl)

#xx<-list.files("/data/sdm_rbq/temp",full=TRUE)
#o<-rev(order(file.info(xx)$mtime))
#head(xx[o])

####################################################################
###  shrink and or quantize images #################################
####################################################################

### web maps
cl<-makeCluster(20)
registerDoParallel(cl)
"/data/sdm_rbq/temp"
lsp<-list.files("/data/sdm_rbq/temp",full=TRUE,pattern=".png")
lsp<-lsp[-grep("_small.",lsp)]
#lsp<-lsp[grep("Spizella_pallida",lsp)]
foreach(i=lsp,.packages=c("magick")) %dopar% {
  image_read(i) |>
  image_scale("500") |>
  image_quantize(max=100,dither=TRUE) |>
  image_write(gsub(".png","_small.png",i))
}
stopCluster(cl)


### marginal effects
cl<-makeCluster(10)
registerDoParallel(cl)
lsp<-list.files("/data/sdm_rbq/marginaleffects",full=TRUE,pattern=".png")
lsp<-lsp[-grep("_small.",lsp)]
#lsp<-lsp[grep("Dendrocygna",lsp)]
foreach(i=lsp,.packages=c("magick")) %dopar% {
  image_read(i) |>
  image_quantize(max=50,dither=TRUE) |>
  image_write(gsub(".png","_small.png",i))
}
stopCluster(cl)


### pixelation and gif
cl<-makeCluster(10)
registerDoParallel(cl)
lsp<-list.files("/data/sdm_rbq/figures",full=TRUE,pattern=".png")
lsp<-lsp[-grep("_small.",lsp)]
#lsp<-lsp[grep("Dendrocygna",lsp)]
foreach(i=lsp,.packages=c("magick")) %dopar% {
  image_read(i) |>
  image_scale("500") |>
  image_quantize(max=100,dither=TRUE) |>
  image_write(gsub(".png","_small.png",i))
}
stopCluster(cl)


### stab
cl<-makeCluster(20)
registerDoParallel(cl)
lsp<-list.files("/data/sdm_rbq/stab",full=TRUE,pattern=".png")
lsp<-lsp[-grep("_small.",lsp)]
foreach(i=lsp,.packages=c("magick")) %dopar% {
  image_read(i) |>
  image_scale("500") |>
  image_quantize(max=50,dither=TRUE) |>
  image_write(gsub(".png","_small.png",i))
}
stopCluster(cl)


### over
cl<-makeCluster(20)
registerDoParallel(cl)
lsp<-list.files("/data/sdm_rbq/overlap",full=TRUE,pattern=".png")
lsp<-lsp[-grep("_small.",lsp)]
foreach(i=lsp,.packages=c("magick")) %dopar% {
  image_read(i) |>
  image_scale("500") |>
  image_quantize(max=50,dither=TRUE) |>
  image_write(gsub(".png","_small.png",i))
}
stopCluster(cl)
