

library(exiftoolr)
library(stringi)
library(ebirdst)
library(magick)

Sys.setlocale("LC_ALL","English")

if(FALSE){
  
  system("powershell -command \"scp -p rouf1703@pose.vhost33:'/data/sdm_rbq/pics/pics.zip' C:/Users/God/Downloads/images\"",intern=TRUE)
  system("powershell -command \"scp -p rouf1703@pose.vhost33:'/data/sdm_rbq/temp/temp.zip' C:/Users/God/Downloads/images\"",intern=TRUE)
  system("powershell -command \"scp -p rouf1703@pose.vhost33:'/data/sdm_rbq/marginaleffects/marginaleffects.zip' C:/Users/God/Downloads/images\"",intern=TRUE)
  system("powershell -command \"scp -p rouf1703@pose.vhost33:'/data/sdm_rbq/overlap/overlap.zip' C:/Users/God/Downloads/images\"",intern=TRUE)
  system("powershell -command \"scp -p rouf1703@pose.vhost33:'/data/sdm_rbq/stab/stab.zip' C:/Users/God/Downloads/images\"",intern=TRUE)
  system("powershell -command \"scp -p rouf1703@pose.vhost33:'/data/sdm_rbq/figures/figures.zip' C:/Users/God/Downloads/images\"",intern=TRUE)
  system("powershell -command \"scp -p rouf1703@pose.vhost33:'/data/sdm_rbq/graphics/mapSpeciesres.csv' C:/Users/God/Downloads/images\"",intern=TRUE)
 
  system("powershell -command \"Get-ChildItem 'C:/Users/God/Downloads/images' -Filter *.zip | Expand-Archive -DestinationPath 'C:/Users/God/Downloads/images' -Force\"",intern=TRUE)
  

  
  
  png("C:/Users/God/Downloads/rejected.png",width=500,height=500,units="px")
  par(mar=c(0,0,0,0))
  plot(0.5,0.5,axes=FALSE,xaxt="n",yaxt="n",pch=1,lwd=20,cex=80,xlim=0:1,ylim=0:1)
  lines(c(0.17,0.83),c(0.17,0.83),lwd=20)
  dev.off()
  im<-image_read("C:/Users/God/Downloads/rejected.png")
  im<-image_trim(im)
  im<-image_scale(im,"500")
  im<-image_quantize(im,2)
  table(image_raster(im)[,3])
  im<-image_transparent(im,"#ffffffff",fuzz=1)
  rejected<-image_fill(im,"#00000033","+250+1")
  #image_write(im,"C:/Users/God/Downloads/rejected.png")
  #file.show("C:/Users/God/Downloads/rejected.png")
  #sdm<-image_read("C:/Users/God/Downloads/images/Setophaga_petechia_sdm_small.png")
  #plot(sdm)
  #plot(im,add=TRUE)
  #im<-image_composite(sdm,im)
  #image_write(im,"C:/Users/God/Downloads/rejected.png")
  #file.show("C:/Users/God/Downloads/rejected.png")
  
  
  df<-read.csv("C:/Users/God/Downloads/images/mapSpeciesres.csv")
  w<-which(df$reach<0.85 | is.na(df$reach))
  invisible(lapply(w,function(i){
    path<-paste0("C:/Users/God/Downloads/images/",gsub(" ","_",df$species[i]),"_sdm_small.png")  
    if(file.exists(path)){
      sdm<-image_read(path)
      #plot(sdm)
      im<-image_composite(sdm,rejected)
      #plot(im)
      image_write(im,path)
    }
  }))
  
}




#plot(1,1)
#tex<-strwrap("This maps show the relative intensity predicted by the model. Intuitively, the relative intensity can be understood as a relative density. It is a measure of relative abundance.")
#text<-paste(tex,collapse="\n")
#text(par("usr")[1],par("usr")[3],text,adj=c(0,0))

mapSpecies_panel<-"This maps shows the relative intensity predicted by the mapSpecies model. Intuitively, the relative intensity can be understood as a relative density surface where high densities indicate a higher number of observations. It is a measure of relative abundance."
eBird_panel<-"This maps shows the mean number of individual birds detected in an eBird checklist for a period of an hour. It is also a measure of relative abundance. To facilitate visualization, a square root transformation has been applied to the color scale. Follow the link over the map to go to the corresponding eBird species page."
overlap_panel<-"This histogram shows the distribution of overlap values between the mapSpecies and the eBird models. A value of 0 indicates completely different distributions and a value of 1 indicates identical distributions. The histogram shows the distribution of values for each species modeled that reached at least 85% of the stabilization hull. The vertical bar shows the overlap value for the focus species."
observations_panel<-"iNaturalist research-grade observations used for modeling the species distributions."
stabilization_panel<-"Shows the stabilization of surface areas of convex hulls around observations when randomly sampling an increasing number of observations. When the curve gets close to the asymptote, it suggests that observations likely cover the range of species."
predictors_panel<-"Shows the relative intensity predicted when the spatial component is removed from the predictions. This represents the effects when only the predictors are considered (whitin the context of the spatial model)."
spatial_panel<-"Shows the spatial effect of the model. Red shows areas where intensities are higher than expected, conditional on the predictors in the model). Blue areas are where intensities are lower than expected. This value is added to the effect of the predictors to form the linear predictor and than exponentiated to get the predicted intensities."
sd_panel<-"Shows the uncertainty (standard deviation) around predictions on the scale of the linear predictor."
marginal_effects_panel<-"The first graphs show the effect of predictors on the relative abundance. The dark line show the predicted mean intensity and the gray area behind the curve shows the 95% credibe intervals. The last two graphs show the posetior"
animation_panel<-""
pixelation_panel<-""

data(ebirdst_runs)
ebird<-as.data.frame(ebirdst_runs)

species_list<-sapply(strsplit(list.files("C:/Users/God/Downloads/images"),"_|\\."),function(i){
  paste(i[1:2],collapse="_") 
}) |> table()
species_list<-species_list[rev(order(species_list))]
species_list<-sort(names(species_list)[species_list>=8])

lf<-list.files("C:/Users/God/Downloads/images",full=TRUE)
m<-match(species_list,gsub("_pic.png","",basename(lf)))
species_list<-species_list[!is.na(m)]
m<-m[!is.na(m)]
url<-exif_read(lf[m],tags="ImageDescription")$ImageDescription
copyright<-exif_read(lf[m],tags="Copyright")$Copyright
copyright <- iconv(copyright, from="UTF-8", to="cp1252") # https://stackoverflow.com/questions/63469654/decoding-cyrillic-string-in-r
Encoding(copyright) <- "UTF-8"
#copyright<-enc2utf8(copyright)
#copyright<-rlang::chr_unserialise_unicode(enc2utf8(enc2native(copyright)))
#copyright<-gsub("\\) ",")<br>",copyright)
copyright<-gsub(", ","<br>",copyright) 
copyright<-gsub(" \\(","<br>(",copyright)
#copyright<-paste0(url,"<br>",copyright)

species_list<-data.frame(sp=species_list,url=url,copyright=copyright)

m<-match(gsub("_"," ",species_list$sp),ebird$scientific_name)
species_code<-ebird$species_code[m]
species_list$ebirdurl<-paste0("https://science.ebird.org/en/status-and-trends/species/",species_code,"/abundance-map?season=breeding&static=true")
species_list$common_name<-ebird$common_name[m]
start<-format(as.Date(ebird$breeding_start[m]),"%b %d")
end<-format(as.Date(ebird$breeding_end[m]),"%b %d")
species_list$period<-ifelse(is.na(start),"Resident",paste(start,end,sep=" / "))
df<-read.csv("C:/Users/God/Downloads/images/mapSpeciesres.csv")
m<-match(gsub("_"," ",species_list$sp),df$species)
species_list$n<-df$n[m]
species_list$family<-df$fname[m]
species_list<-species_list[order(match(species_list$common,ebird$common)),]
species_list<-species_list[order(factor(species_list$family,levels=unique(species_list$family))),] # some in ebird appear unordered

#species_list$family<-sort(rep(LETTERS[1:20],length.out=nrow(species_list)))
#species_list<-species_list[sample(1:5,20,replace=TRUE),]
#species_list<-species_list[order(species_list$family),]
toc_species<-rep(unique(species_list$family),times=rle(species_list$family)[[1]]+1)
toc_ref<-rep(unique(species_list$sp),times=as.integer(!duplicated(species_list$family))+1)
toc_common<-rep(unique(species_list$common),times=as.integer(!duplicated(species_list$family))+1)
toc_species<-ifelse(duplicated(toc_ref,fromLast=TRUE),paste0("<b>",toc_species,"</b>"),paste0("&nbsp;&nbsp;",toc_common))
head(data.frame(toc_ref,toc_species))



#src<-"https://res.cloudinary.com/dphvzalf9/image/upload"
src<-"images"

species<-function(sp,url,copyright,ebirdurl,common,period,n){
  
  paste0("

<hr class=\"vspace\"> 
<h2 id=\"",sp,"\" class=\"h2\">",common,"</h2>  
<!-- <section class=\"section\"></section> -->
<div class=\"subheader\"></div>


<section class=\"section\">
  <header>
    <div class=\"top-left\"></div>
    <div class=\"top-left\">mapSpecies</div>
    <div class=\"top-left\">
      <a target=\"_blank\" href=\"",ebirdurl,"\">eBird &#8599</a>
    </div>
    <div class=\"top-left\">Overlap (I)</div>
  </header>
  <div class=\"rowshow\">
    <div class=\"col3\" style=\"width: 20vw; border: 0px solid red;\">
      <!-- <a target=\"_blank\" href=\"",url,"\"> -->
      <div class=\"container\" style=\"border: 0px solid red;\">
        <a class=\"aim\" target=\"_blank\" href=\"",url,"\">
          <div class=\"infotop\">",period,"<br>n = ",n,"</div>
          <!-- <a class=\"aim\" target=\"_blank\" href=\"",url,"\"> -->
          <img style=\"height: 12vw; padding-top: 5%; padding-left: 0%; border: 0px solid green;\" src=\"",file.path(src,sp),"_pic.png\" alt=\"\">
        <div class=\"middle\" style=\"border: 0px solid red;\">
          <div class=\"text\" style=\"border: 0px solid green;\">",copyright,"</div>
        </div>
        <!-- </a> -->
        <div class=\"infobottom\"><i>",gsub("_"," ",sp),"</i></div>
        </a>
      </div>
    </div>
    <div class=\"col4\">
      <img style=\"height: 20vw; padding: 0px;\" src=\"",file.path(src,sp),"_sdm_small.png\" alt=\"\">
    </div>
    <div class=\"col3\">
      <img style=\"height: 20vw; padding: 0px;\" src=\"",file.path(src,sp),"_ebird_small.png\" alt=\"\">
    </div>
    <div class=\"col4\">
      <img style=\"height: 16vw; padding-top: 20%;\" src=\"",file.path(src,sp),"_over.png\" alt=\"\">
    </div>
  </div>
</section>





<section class=\"section\">
  <header>
    <button class=\"showmore\" onclick=\"showmorefirst('",sp,"','panel1","')\">&nbsp&nbsp&nbsp&nbspObservations and stabilization &#8628</button>
    <button class=\"showmore\" onclick=\"showmorefirst('",sp,"','panel2","')\">Predictors only and spatial effect &#8628</button>
  </header>
  <div class=\"row\" id=\"",sp,"panelfirst\">
    <div class=\"col1\" id=\"",sp,"panel1\">
      <img style=\"height: 17vw; padding: 0px;\" src=\"",file.path(src,sp),"_obs_small.png\" alt=\"\">
      <img style=\"height: 17vw; padding: 0px;\" src=\"",file.path(src,sp),"_stab_small.png\" alt=\"\">
    </div>
    <div class=\"col2\" id=\"",sp,"panel2\">
      <img style=\"height: 17vw; padding: 0px;\" src=\"",file.path(src,sp),"_pred_small.png\" alt=\"\">
      <img style=\"height: 17vw; padding: 0px;\" src=\"",file.path(src,sp),"_spatial_small.png\" alt=\"\">
    </div>
  </div>
</section>


<section class=\"section\">
  <header>
    <button class=\"showmore\" onclick=\"showmoresecond('",sp,"','panel3","')\">&nbsp&nbsp&nbsp&nbspSD, pixelation and animated uncertainty &#8628</button>
    <button class=\"showmore\" onclick=\"showmoresecond('",sp,"','panel4","')\">Marginal effects of predictors &#8628</button>
  </header>
  <div class=\"row\" id=\"",sp,"panelsecond\">
    <div class=\"col1\" id=\"",sp,"panel3\">
      <img style=\"height: 17vw; padding: 0px;\" src=\"",file.path(src,sp),"_sd_small.png\" alt=\"\">
      <img style=\"height: 17vw; padding: 0px;\" src=\"",file.path(src,sp),"_pixelation_small.png\" alt=\"\">
      <img style=\"height: 17vw; padding: 0px; opacity: 0;\" src=\"",file.path(src,sp),"_sd_small.png\" alt=\"\">
      <img style=\"height: 17vw; padding: 0px;\" src=\"",file.path(src,sp),"_gif.gif\" alt=\"\">
    </div>
    <div class=\"col2\" id=\"",sp,"panel4\">
      <img style=\"height: 29vw; padding: 0px;\" src=\"",file.path(src,sp),"_me_small.png\" alt=\"\">
    </div>
  </div>
</section>


  
  ")
}


css<-function(){cat(paste0("

<!DOCTYPE html>
<html>
<head>
<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">
<link href='https://fonts.googleapis.com/css2?family=Roboto+Mono:wght@400' rel='stylesheet'>

<style>

html {
  /* max-width: 1950px !important; */
  background-color: white;
}

body {
  background-color: white;
  max-width: 95vw;
  margin-left: auto;
  margin-right: auto;
  font-family: 'Roboto Mono';
  font-size: 2.25vmin;
}

div.left {
  position: -webkit-sticky;
  position: sticky;
  top: 3vh;
  bottom: 3vh;
  width: 18%; 
  height: 94vh; 
  float: left; 
  display: inline-block; 
  padding-right: 0.0%; 
  overflow-y: scroll;
}



div.right {
  width: 80%; 
  float: left; 
  display: inline-block; 
  padding-left: 2%;
}


a:link, a:visited, a:hover, a:active {
 color: #000000;
 text-decoration: none;
}

a:hover{
  opacity: 0.50;
}

a.toc {
  font-size: 1.75vmin;
  color: black;
}

.aim:hover{
  opacity: 1;
}

a.text{
  font-size: 2.25vmin;
  font-weight: 600;
}

a.text:link, a.text:visited, a.text:hover, a.text:active {
 color: #000000;
 text-decoration: none;
}

a.text:hover{
  opacity: 0.50;
}

.vspace {
  margin-top: 2vh;
  margin-bottom: 12vh;
  padding:0;
  border: 0;
  background-color: white;
}
.subheader {
  height: 1vh;
  margin: 0vh;
  padding: 0vh;
  border-left: 0.5vmin solid seagreen;
}

.infotop {
  padding-top: 10%;
  font-size: 1.75vmin;
}

.infobottom {
  padding-top: 5%;
  font-size: 1.75vmin;
}

.container {
  position: relative;
  text-align: center;
  color: black;
}

.container2 {
  position: relative;
  text-align: center;
  color: black;
}

.top-left {
  display: flex;
  flex: 1;
  padding: 0vmin;
  margin: 0vmin;
  /* border: 1px solid red; */
}

h1 {
  color: seagreen;
}

.h2 {
  margin: 0vh;
  padding: 0vh;
  color: seagreen;
  /* border-bottom: 2px solid seagreen; */
  /* border-right: 2px solid seagreen; */
}

.section{
  border-left: 0.5vmin solid seagreen;
  /* border-bottom: 2px solid seagreen; */
  /* border-right: 2px solid seagreen; */
}

/* .hide {
  display: none;
  padding-top: 2vw;
  padding-bottom: 0vw;
} 

.show {
  padding-top: 0vw;
  font-weight: 1200;
  font-size: 2.5vmin;
}

.show:hover + .hide {
  display: block;
} */

.shown {
  display: block;
  padding-top: 0vw;
  border-left: 0.5vmin solid seagreen;
  <!-- border-right: 2px solid seagreen; -->
}

header, .row {
  display: flex;  /* aligns all child elements (flex items) in a row */
}

header {
  display: flex;  /* aligns all child elements (flex items) in a row */
  padding: 0vmin;
  margin: 0vmin;
}

.row {
  display: none;  /* aligns all child elements (flex items) in a row */
  padding-top: 1vh;
  padding-bottom: 1vh;
}

.rowshow {
  display: flex;  /* aligns all child elements (flex items) in a row */
  padding-top: 0vh;
  padding-bottom: 1vh;
}

.showmore {
  background: none;
  border: none; /* none */
  padding-bottom: 3vh;
  text-align: center;
  text-decoration: none;
  display: flex;
  flex: 1;
  font-size: 2.25vmin;
  font-weight: 1200;
  font-family:'Roboto Mono'; 
  cursor: pointer;
  height: 2vmin;
}

.showmore:hover {
  opacity: 0.50;
  filter: alpha(opacity=100);
}

.col1, .col2 {
  display: flex;
  flex: 1;
  visibility: hidden;
  flex-flow: wrap;
}

.col3 .col4 {
  display: flex;
  flex: 1;
  visibility: show;
  transition: transform .2s;
}

/* .col4:hover {
  -ms-transform: scale(1.1); /* IE 9 */
  -webkit-transform: scale(1.1); /* Safari 3-8 */
  transform: scale(1.1); 
} */



.ID {
  display: none;
}

#####
.container {
  position: relative;
  width: 50%;
  border: 0px solid red;
}

.container2 {
  position: relative;
  width: 100%;
  border: 0px solid red;
  opacity: 1;
}

.container3 {
  position: relative;
  /* width: 50%; */
  /* border: 1px solid red; */
}

.image {
  opacity: 1;
  display: block;
  width: 100%;
  height: auto;
  transition: .2s ease;
  backface-visibility: hidden;
}

.middle {
  transition: .2s ease;
  opacity: 0;
  position: absolute;
  top: 70%;
  left: 50%;
  width: 60%;
  transform: translate(-50%, -50%);
  -ms-transform: translate(-50%, -50%);
  text-align: center;
}

.middle2 {
  opacity: 1;
  position: absolute;
  background-color: #FFFFFF00;
  top: 25%;
  left: 15%;
  width: 80%;
  text-align: left;
}

.middle3 {
  opacity: 1;
  position: absolute;
  background-color: #FFFFFF00;
  top: 17%;
  left: 2%;
  width: 98%;
  text-align: left;
  /* border: 1px solid red; */
}

.container:hover .image {
  opacity: 0.3;
}

.container:hover .middle {
  opacity: 1;
}

.text {
  background-color: #FFFFFFAA;
  color: #000000;
  font-size: 1.25vmin;
  padding: 0.25vmin;
  border-radius: 1vmin;
  line-height: 1;
}

.text2 {
  background-color: #FFFFFF00;
  color: #00000033;
  opacity: 1;
  font-size: 1.75vmin;
  padding: 1vmin;
  border-radius: 1vmin;
  line-height: 1.1;
  /* text-align:justify; */
  text-justify:inter-character;
}

.container2:hover .image2 {
  opacity: 0.4;
}

.container2:hover .image3 {
  opacity: 0.2;
}

.container2:hover .text2 {
  color: #000000;
  background-color: #FFFFFF55;
}

.container3:hover .image2 {
  opacity: 0.4;
}

.container3:hover .image3 {
  opacity: 0.2;
}

.container3:hover .text2 {
  color: #000000;
  background-color: #FFFFFF55;
}

</style>
</head>
<body>
<div style=\"display:inline-block; width:100%;\">
<div class=\"left\">"
,paste(paste0("<a class=\"toc\" href=\"#",toc_ref,"\">",toc_species,"</a>"),collapse="<br>\n"),
"</div>
<div class=\"right\">
<h1>SDM from mapSpecies and comparison with eBird</h1>

<p>This page displays the species distribution models obtained with mapSpecies and associated results along with the abundance models from <a class=\"text\" target=\"blank\" href=\"https://ebird.org\">eBird</a> for comparison. All observations used for the mapSpecies models come from <a class=\"text\" target=\"blank\" href=\"https://www.inaturalist.org\">iNaturalist</a> and were obtained through the <a class=\"text\" target=\"blank\" href=\"https://doi.org/10.15468/ab3s5x\">iNaturalist Research-grade Observations</a> dataset hosted on <a class=\"text\" target=\"blank\" href=\"https://www.gbif.org\">GBIF</a>. The eBird models are part of <a class=\"text\"  target=\"blank\" href=\"https://science.ebird.org/en/status-and-trends\">The Status and Trends product</a> <a class=\"text\" target=\"blank\" href=\"https://doi.org/10.2173/ebirdst.2021\">(Fink et al. 2022)</a> and can be accessed through the <a class=\"text\" target=\"blank\" href=\"https://cornelllabofornithology.github.io/ebirdst/index.html\">ebirdst</a> R package. The abundance models from eBird are available online and a link to the ebird model is provided in each species account. The models are either for the breeding period or for the whole year when a species is a resident. The breeding period used for the models is the same as the one used by eBird.</p>

<hr class=\"vspace\"> 
<h2 id=\"Results\" class=\"h2\">Paper figures</h2> 
<figure>
  <img style=\"height: 75vh; padding: 0px;\" src=\"",file.path(src,"visual_abstract.png"),"\" alt=\"\">
  <figcaption>Figure 1.</figcaption>
</figure><br>
<figure>
  <img style=\"height: 50vh; padding: 0px;\" src=\"",file.path(src,"overlap.png"),"\" alt=\"\">
  <figcaption>Figure 2.</figcaption>
</figure><br>
<figure>
  <img style=\"height: 45vh; padding: 0px;\" src=\"",file.path(src,"spatial_effect.png"),"\" alt=\"\">
  <img style=\"height: 45vh; padding: 0px;\" src=\"",file.path(src,"effort_effect.png"),"\" alt=\"\">
  <figcaption>Figure 3.</figcaption>
</figure><br>
<figure>
  <img style=\"height: 100vh; padding: 0px;\" src=\"",file.path(src,"bird_maps.png"),"\" alt=\"\">
  <figcaption>Figure 5.</figcaption>
</figure><br>
<figure>
  <img style=\"height: 100vh; padding: 0px;\" src=\"",file.path(src,"results_detailed.png"),"\" alt=\"\">
  <figcaption>Figure 6.</figcaption>
</figure><br>


<hr class=\"vspace\"> 
<h2 id=\"Explanations\" class=\"h2\">See below and hover over figures for quick explanations</h2>  
<!-- <section class=\"section\"></section> -->
<div class=\"subheader\"></div>


<section class=\"section\">
  <header>
    <div class=\"top-left\"></div>
    <div class=\"top-left\">mapSpecies</div>
    <div class=\"top-left\">
      <a target=\"_blank\" href=\"https://science.ebird.org/en/status-and-trends/species/bbwduc/abundance-map?season=breeding&static=true\">eBird &#8599</a>
    </div>
    <div class=\"top-left\">Overlap (I)</div>
  </header>
  <div class=\"rowshow\">
    <div class=\"col3\" style=\"width: 20vw; border: 0px solid red;\">
      <!-- <a target=\"_blank\" href=\"https://www.inaturalist.org/observations/142367503\"> -->
      <div class=\"container\" style=\"border: 0px solid red;\">
        <a class=\"aim\" target=\"_blank\" href=\"https://www.inaturalist.org/observations/142367503\">
          <div class=\"infotop\">Breeding period<br>No. of observations used</div>
          <!-- <a class=\"aim\" target=\"_blank\" href=\"https://www.inaturalist.org/observations/142367503\"> -->
          <img style=\"height: 12vw; padding-top: 5%; padding-left: 0%; border: 0px solid green;\" src=\"",file.path(src,"Dendrocygna_autumnalis"),"_pic.png\" alt=\"\">
        <div class=\"middle\" style=\"border: 0px solid red;\">
          <div class=\"text\" style=\"border: 0px solid green;\">(c) Author of photo<br>Licence<br>Click to see observation<br>on iNaturalist</div>
        </div>
        <!-- </a> -->
        <div class=\"infobottom\"><i>Scientific name</i></div>
        </a>
      </div>
    </div>
    <div class=\"col4\">
    <div class=\"container2\">
      <img class=\"image2\" style=\"height: 20vw; padding: 0px;\" src=\"",file.path(src,"Dendrocygna_autumnalis"),"_sdm_small.png\" alt=\"\">
      <div class=\"middle2\" style=\"border: 0px solid red;\">
      <div class=\"text2\" style=\"border: 0px solid green;\">",mapSpecies_panel,"</div>
      </div>
    </div>
    </div>
    <div class=\"col3\">
    <div class=\"container2\">
      <img class=\"image2\" style=\"height: 20vw; padding: 0px;\" src=\"",file.path(src,"Dendrocygna_autumnalis"),"_ebird_small.png\" alt=\"\">
      <div class=\"middle2\" style=\"border: 0px solid red;\">
      <div class=\"text2\" style=\"border: 0px solid green;\">",eBird_panel,"</div>
      </div>
    </div>
    </div>
    <div class=\"col4\">
    <div class=\"container2\">
      <img class=\"image3\" style=\"height: 16vw; padding-top: 20%;\" src=\"",file.path(src,"Dendrocygna_autumnalis"),"_over.png\" alt=\"\">
      <div class=\"middle2\" style=\"border: 0px solid red;\">
      <div class=\"text2\" style=\"border: 0px solid green;\">",overlap_panel,"</div>
      </div>
    </div>
    </div>
  </div>
</section>





<section class=\"section\">
  <header>
    <button class=\"showmore\" onclick=\"showmore('Explanations','panel1')\">&nbsp&nbsp&nbsp&nbspObservations and stabilization &#8628</button>
    <button class=\"showmore\" onclick=\"showmore('Explanations','panel2')\">Predictors only, spatial effect and SD &#8628</button>
  </header>
  <div class=\"row\" id=\"Explanationspanel\" style=\"display: flex;\">
    <div class=\"col1\" id=\"Explanationspanel1\" style=\"visibility: visible; display: flex;\">
      <div class=\"container3\">
      <img class=\"image2\" style=\"height: 13vw; padding: 0px;\" src=\"",file.path(src,"Dendrocygna_autumnalis"),"_obs_small.png\" alt=\"\">
      <div class=\"middle3\" style=\"border: 0px solid red;\">
      <div class=\"text2\" style=\"border: 0px solid green;\">",observations_panel,"</div>
      </div>
      </div>
      <div class=\"container3\">
      <img class=\"image3\" style=\"height: 13vw; padding: 0px;\" src=\"",file.path(src,"Dendrocygna_autumnalis"),"_stab.png\" alt=\"\">
      <div class=\"middle3\" style=\"border: 0px solid red;\">
      <div class=\"text2\" style=\"border: 0px solid green;\">",stabilization_panel,"</div>
      </div>
      </div>
    </div>
    <div class=\"col2\" id=\"Explanationspanel2\" style=\"visibility: visible; display: flex;\">
      <div class=\"container3\">
      <img class=\"image2\" style=\"height: 13vw; padding: 0px;\" src=\"",file.path(src,"Dendrocygna_autumnalis"),"_pred_small.png\" alt=\"\">
      <div class=\"middle3\" style=\"border: 0px solid red;\">
      <div class=\"text2\" style=\"border: 0px solid green;\">",predictors_panel,"</div>
      </div>
      </div>
      <div class=\"container3\">
      <img class=\"image2\" style=\"height: 13vw; padding: 0px;\" src=\"",file.path(src,"Dendrocygna_autumnalis"),"_spatial_small.png\" alt=\"\">
      <div class=\"middle3\" style=\"border: 0px solid red;\">
      <div class=\"text2\" style=\"border: 0px solid green;\">",spatial_panel,"</div>
      </div>
      </div>
      <div class=\"container3\">
      <img class=\"image2\" style=\"height: 13vw; padding: 0px;\" src=\"",file.path(src,"Dendrocygna_autumnalis"),"_sd_small.png\" alt=\"\">
      <div class=\"middle3\" style=\"border: 0px solid red;\">
      <div class=\"text2\" style=\"border: 0px solid green;\">",sd_panel,"</div>
      </div>
      </div>
    </div>
  </div>
</section>

<hr class=\"vspace\"> 



"))}


script<-function(){cat(paste0("

<script>
    var buttons = document.getElementsByClassName('showmore');
    var divs = document.getElementsByClassName('ID');
    
    for (var i = 0; i < buttons.length; i++) {
      var butt = buttons[i];
      var div = divs[i];
      // and attach our click listener for this image.
      buttsssssss.onclick = function(evt) {
        console.log(evt);
        if (div.style.display == \"flex\") {
          div.style.display = \"none\";
        } else {
          div.style.display = \"flex\";
        }
      }
    }
    
    
    function showmorefirst(sp,panel) {
      var id = sp+panel
      var x = document.getElementById(id);
      var div = document.getElementById(sp+\"panelfirst\");
      var div1 = document.getElementById(sp+\"panel1\");
      var div2 = document.getElementById(sp+\"panel2\");
      var status = \"hidden\";
      if (x.style.visibility == \"visible\") {
        x.style.visibility = \"hidden\";
      } else {
        x.style.visibility = \"visible\";
        status = \"visible\";
      }
      if (div1.style.visibility == \"visible\" || div2.style.visibility == \"visible\" || status == \"visible\") {
        div.style.display = \"flex\";
      } else {
        div.style.display = \"none\";
      }
    }
    
    function showmoresecond(sp,panel) {
      var id = sp+panel
      var x = document.getElementById(id);
      var div = document.getElementById(sp+\"panelsecond\");
      var div3 = document.getElementById(sp+\"panel3\");
      var div4 = document.getElementById(sp+\"panel4\");
      var status = \"hidden\";
      if (x.style.visibility == \"visible\") {
        x.style.visibility = \"hidden\";
      } else {
        x.style.visibility = \"visible\";
        status = \"visible\";
      }
      if (div3.style.visibility == \"visible\" || div4.style.visibility == \"visible\" || status == \"visible\") {
        div.style.display = \"flex\";
      } else {
        div.style.display = \"none\";
      }
    }
    
</script>

"))}
  




con <- file("C:/Users/God/Downloads/website.html", open = "w+b")#, encoding = "UTF-8")
sink(con)


css()  
invisible(lapply(1:nrow(species_list),function(i){
  ans<-species(species_list$sp[i],species_list$url[i],species_list$copyright[i],species_list$ebirdurl[i],species_list$common_name[i],species_list$period[i],species_list$n[i])
  stri_write_lines(ans,con=con)
}))
script()
cat(paste0("
</body>
</html>
"))


close(con)

file.show("C:/Users/God/Downloads/website.html")

if(FALSE){
 
  
  getCC0links<-function(species,license=c("cc0")){  
    sp<-gsub(" ","%20",species)
    cc<-paste(license,collapse="0%2C")
    urlsearch<-paste0("https://api.inaturalist.org/v1/taxa?q=",sp,"&order=desc&order_by=observations_count")
    x<-fromJSON(urlsearch)$results
    taxonid<-x$id[1]
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
    organim<-image_trim(organism,fuzz=0)
    organism<-image_fill(organism,"none",point="+5+5",fuzz=2)
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
  
  library(jsonlite)
  library(magick)
  
  path<-"C:/Users/God/Downloads/images"
  
  #df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
  #sps<-c("Buteo jamaicensis","Buteo lagopus","Buteo platypterus")
  #sps<-sort(df$species)#[1:18]
  #sps<-sps[165:length(sps)]
  i<-c("Acanthis hornemanni")
  sps<-i#c("Columba livia")
  
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
  
  
  
  con <- file("C:/Users/God/Downloads/test.html",open="w+b")#, open = "wt", encoding = "latin1")
  sink(con)
  
  css() 
  ans<-species(species_list$sp[1],species_list$url[1],species_list$copyright[1])
  stri_write_lines(ans,con=con)
  script()
  cat(paste0("
</div>
</div>
</body>
</html>
"))
  
  close(con)
  
  file.show("C:/Users/God/Downloads/test.html")
  
  
  
  
}



