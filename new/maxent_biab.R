
# https://stackoverflow.com/questions/34624002/r-error-java-lang-outofmemoryerror-java-heap-space
options(java.parameters = "-Xmx40g")

library(data.table)
library(mapSpecies)
library(sf)
library(terra)
library(concaveman)
library(FNN)
library(predicts)
library(rmapshaper)
library(s2)
library(ewlgcpSDM)
library(FRutils)
library(future)
library(future.apply)

options(vsc.dev.args = list(width = 1000, height = 800))


flex_buffer<-function(
    obs,
    nline=48,
    dist=c(500,1500)
){
  # https://r-spatial.github.io/sf/articles/sf7.html#buffers-1

  radius<-s2_earth_radius_meters()*pi*2*0.25*((90-nline)/90)
  npole<-s2_buffer_cells(as_s2_geography("POINT(-100 90)"),distance=radius,max_cells=10000) # visible half
  g<-as_s2_geography(TRUE)
  north<-s2_intersection(npole, g)
  north<-st_transform(st_as_sfc(north),st_crs(na))



  buff1<-st_union(st_buffer(obs[lengths(st_intersects(obs,north))==0L,],dist[1]))
  buff2<-st_union(st_buffer(obs[lengths(st_intersects(obs,north))>0L,],dist[2]))

  bnorth<-st_cast(st_intersection(buff2,north),"POLYGON")
  o<-st_intersects(bnorth,na)
  bnorth<-bnorth[lengths(o)>0]
  #bsouth<-st_cast(st_difference(buff1,north),"POLYGON")
  bsouth<-st_cast(buff1,"POLYGON")
  o<-st_intersects(bsouth,na)
  bsouth<-bsouth[lengths(o)>0]
  buff<-st_union(bnorth,bsouth) |> st_union()
  #plot(st_geometry(na))
  #plot(st_geometry(buff),col=adjustcolor("black",0.15),border=NA,add=TRUE)
  #plot(st_geometry(obs),add=TRUE)
  #plot(st_geometry(bnorth),col="blue",add=TRUE)
  #plot(st_geometry(bsouth),col="red",add=TRUE)
  #lpols<-lapply(buff, function(x) x[1])
  #buff<-st_multipolygon(lpols) # remove holes https://github.com/r-spatial/sf/issues/609
  buff
}



colmean<-colo.scale(1:200,c("#CCCCCC","#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00","darkorange"))
prj<-"+proj=lcc +lat_0=47 +lon_0=-75 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

predictors<-rast("/data/predictors_sdm/predictors.tif")
predictors<-predictors[[c("tmean","tmax","trange","prec","elevation","truggedness","longitude","latitude",grep("_esa",names(predictors),value=TRUE))]]
predictors<-aggregate(predictors,4,na.rm=TRUE)

#scale(predictors)

na<-st_read(list.files("/data/predictors_sdm",pattern="na.shp",full=TRUE))
na<-st_transform(na,st_crs(predictors))
coast<-st_union(na)
nalakes<-st_read("/data/predictors_sdm/nalakes.gpkg")
nalakes<-st_transform(nalakes,st_crs(na))
areas<-as.numeric(st_area(nalakes))
ke<-paste(c("Manicouagan","-Jean"),collapse="|")
nalakes<-nalakes[unique(c(rev(order(areas))[1:50],grep(ke,nalakes$name_fr))),]


gbif<-rast("https://object-arbutus.cloud.computecanada.ca/bq-io/io/gbif_heatmaps/gbif_plants_density_06-2022.tif")
gbif<-crop(gbif,vect(st_transform(na,4326)))
gbif<-mask(gbif,vect(st_transform(na,4326)))
gbif1<-subst(gbif,NA,0)
gbif<-rast("https://object-arbutus.cloud.computecanada.ca/bq-io/io/gbif_heatmaps/gbif_reptiles_density_06-2022.tif")
gbif<-crop(gbif,vect(st_transform(na,4326)))
gbif<-mask(gbif,vect(st_transform(na,4326)))
gbif2<-subst(gbif,NA,0)

gbif<-gbif1+gbif2
gbif<-gbif1
#gbif<-aggregate(gbif,5,fun=sum,na.rm=TRUE)


# https://coleo.biodiversite-quebec.ca/apps/tableau-forets/
trees<-rast("https://object-arbutus.cloud.computecanada.ca/bq-io/io/forets-cc-landis/baseline_BudwormBaselineFire_ACER.RUB_0_merged.tif")
trees[is.na(trees)]<-0
trees<-mask(trees,vect(st_transform(na[na$NAME_1=="Québec",],st_crs(trees))))
trees<-aggregate(trees,10,na.rm=TRUE)
trees<-project(trees,crs(na))
trees<-focal(trees,w=5,fun="mean",na.rm=TRUE)
plot(trees)
plot(st_geometry(na),add=TRUE)


#region<-na[na$NAME_1%in%c("Québec"),]
region<-na[na$NAME_1%in%c("Québec","New Brunswick","Maine","Vermont","New Hampshire","New York","Ontario"),]
labrador<-ms_explode(na[na$NAME_1%in%c("Newfoundland and Labrador"),])
labrador<-labrador[which.max(st_area(labrador)),]
region<-rbind(region,labrador)
regionll<-st_transform(region,4326)
wkt<-ms_simplify(region,0.01)
wkt<-st_union(wkt)
wkt<-st_transform(wkt,4326)
wkt<-st_as_text(wkt)

gbifreg<-mask(crop(gbif,vect(regionll)),vect(regionll))
predictorsreg<-mask(crop(predictors,vect(region)),vect(region))

obs<-rgbif::occ_data(scientificName = "Fagus grandifolia", hasCoordinate = TRUE,limit=5000,geometry=wkt)$data
rem<-which(obs$coordinateUncertaintyInMeters>30000)
if(any(rem)){
  obs<-obs[-rem,]
}
obs<-st_as_sf(as.data.frame(obs),coords=c("decimalLongitude","decimalLatitude"),crs=4326)
obs<-st_transform(obs,st_crs(region))

fg<-predictors[[1]]
values(fg)<-1:ncell(fg)
e<-terra::extract(fg,vect(obs))
dups<-duplicated(data.frame(cell=e[,2],observer=obs$recordedBy))
table(dups)
obs<-obs[!dups,]



# Apply filter grid


#obs<-obs[obs$recordedBy!="Louis Imbeau",]

plot(st_geometry(region))
plot(st_geometry(na),col="grey90",add=TRUE)
plot(st_geometry(region),col="white",add=TRUE)
plot(st_geometry(obs),add=TRUE)


buff<-flex_buffer(obs,nline=48,dist=c(500,500))
tbplus<-st_as_sf(st_sample(st_difference(region,buff),10000))
st_geometry(tbplus)<-"geometry"
#plot(st_geometry(tbplus),add=TRUE)
p<-spatSample(gbifreg,500000,as.points=TRUE,method="weights")
tb<-st_as_sf(p) |> st_transform(st_crs(obs))
tb<-tb[,"geometry"]
tb<-rbind(tb,tbplus)
#tb<-st_sample(region,50000)

cloglog<-function(x){
  g<-global(x,"range",na.rm=TRUE)
  y<-log(-log(((g[1,2]-x)/(g[1,2]-g[1,1]))))
  y<-1/(1+exp(-y))
  y/global(y,"max",na.rm=TRUE)[1,1]
}

# fit model
#me <- MaxEnt(predictorsreg, vect(obs),vect(tb),args=c("-q","-l","-p"),removeDuplicates=TRUE,silent=FALSE)#, vect(a))
#g<-names(predictorsreg)
g<-unique(grep("2",vs,value=TRUE,invert=TRUE))
#g<-g[c(1,2,3,6)]
test<-predictorsreg[[g]]
me <- MaxEnt(test, vect(obs),vect(tb),removeDuplicates=FALSE,silent=FALSE)#, vect(a))
mer <- mask(predict(me, test,args=c("outputformat=raw")),vect(na))
plot(mer,mar=c(0.5,0.5,0.5,8))
#plot(cloglog(mer),mar=c(0.5,0.5,0.5,8))
plot(st_geometry(region),add=TRUE)
plot(st_geometry(nalakes),col="white",border=NA,add=TRUE)
plot(st_geometry(na),add=TRUE)
points(obs,pch=1,cex=0.5,col=adjustcolor("black",0.5))
#points(tb,pch=1,cex=0.5,col=adjustcolor("black",0.5))


#hist(params$predictors$builtup_esa)


#############################
#############################
#############################
#############################
#############################
#############################


library(ewlgcpSDM)

domain<-st_sample(st_buffer(region,50),5000)
domain <- inla.nonconvex.hull(as(domain,"Spatial"),convex = -0.015,resolution=75)

pedge<-0.004
edge<-min(c(diff(st_bbox(region)[c(1,3)])*pedge,diff(st_bbox(region)[c(2,4)])*pedge))
edge

mesh <- inla.mesh.2d(loc.domain = NULL,
                     max.edge = c(edge,edge*3),
                     min.angle = 21,
                     cutoff = edge/1,
                     offset = c(edge,edge*3),
                     boundary = domain,#inla.mesh.segment(domain$loc),
                     crs = st_crs(region))


bg<-rbind(obs[,"geometry"],tb)

#plan(multicore,workers=10)
params<-dmesh_mesh(mesh)
params<-dmesh_weights(params,region)
params<-dmesh_predictors(params,predictorsreg)
params<-dmesh_effort(params,obs=obs,background=bg,buffer=buff,adjust=FALSE)
#plan(sequential)

params$effort$nobs<-ifelse(params$predictors$builtup_esa>0.5,0,params$effort$nobs)
params$effort$nbackgroundadjusted<-ifelse(params$predictors$builtup_esa>0.5,0,params$effort$nbackgroundadjusted)

#params$effort$nobs<-as.integer(params$effort$nobs>0)
#hist(params$effort$nobs/params$effort$nbackground)

dm<-cbind(params$dmesh,params$predictors,params$effort)
dm$ratio<-dm$nobs/dm$nbackground
#plot(dm["nobs"],border=NA,reset=FALSE)
#plot(st_geometry(na),border="white",add=TRUE)
#plot(st_geometry(nalakes),col="white",border=NA,add=TRUE)
#dev.off()

plot(st_geometry(dm),border=NA,col=ifelse(dm$builtup_esa>0.5,"red","blue"))
plot(st_geometry(na),border="white",add=TRUE)


#eff<-ifelse((params$effort$nobs/params$effort$nbackground)>0.1,10,0)
#eff[is.na(eff)]<-0
#params$effort$nbackgroundadjusted<-params$effort$nbackgroundadjusted+eff

vs<-c("tmax","tmax2","latitude","latitude2","longitude","longitude2","conifers_esa","deciduous_esa","conifers_esa","crop_esa","mixed_esa","water_esa","truggedness","wettree_esa","grass_esa","wetherbaceous_esa")
#vs<-paste0(sort(rep(grep("2",vs,value=TRUE,invert=TRUE),2)),c("","2"))
mf<-paste(vs,collapse="+")

m<-ewlgcp(
  formula=formula(paste("y~",mf)),
  dmesh=params,
  effort = TRUE,
  adjust = FALSE,
  buffer = TRUE,
  orthogonal = TRUE,
  prior.beta = NULL,
  prior.range = c(50,0.01),
  prior.sigma = c(0.000001,NA),
  smooth = 3/2
)

sdm<-ewlgcpSDM::map(model=m,
         dmesh=params,
         dims=c(1500,1500),
         region = region
     )

sdm<-mask(sdm,vect(na))
crs(sdm)<-crs(na)


plot(mask(sdm$mean,region),plg=list(size=c(0.8,1)))
plot(st_geometry(na),lwd=0.5,add=TRUE)
plot(st_geometry(nalakes),col="white",lwd=0.2,add=TRUE)
plot(st_geometry(obs),add=TRUE)


par(mfrow=c(2,2))
plot(crop(mer,ext(trees)));plot(st_geometry(st_transform(na,st_crs(trees))),add=TRUE)
#plot(mask(project(mer,trees),trees))
plot(trees);plot(st_geometry(na),add=TRUE)
plot(crop(mask(sdm$mean,region),trees),plg=list(size=c(0.8,1)));plot(st_geometry(na),add=TRUE)
plot(st_geometry(nalakes),col="white",lwd=0.2,add=TRUE)
par(mfrow=c(1,1))




##########################
##########################
##########################
##########################
##########################


  dm<-st_transform(dmesh$dmesh,st_crs(predictors))
  cores<-nbrOfWorkers() # get nbr of workers from the chosen plan
  chunks <- split(1:nrow(dm), rep(1:cores, each=ceiling(nrow(dm)/cores))[1:nrow(dm)])
  options(future.globals.maxSize = 1000 * 1024 ^ 2)
  res<-future_lapply(chunks,function(chunksi){
    t(exact_extract(predictors,
                    dm[chunksi,],
                    fun = function(values, coverage_fraction){
                      colSums(as.matrix(values) * coverage_fraction,
                              na.rm = TRUE) / sum(coverage_fraction)
                    },
                    force_df = FALSE,
                    progress = TRUE))
    #res<-t(res)
  })
  #plan(sequential)
  res<-as.data.frame(do.call("rbind",res))

  #dm<-cbind(dm,res)
  #plot(dm["conifers_esa"])
  #plot(predictors$tmean)


  ### Fill mesh where values are missing
  xy<-st_coordinates(st_centroid(dm))
  xy<-cbind(xy,res)
  xynotna<-xy |> as.data.table() |> na.omit() |> as.matrix()
  nn<-knnx.index(xynotna[,1:2],xy[,1:2],k=1)
  res<-as.data.frame(xynotna[nn,-(1:2)])





rivers <- ne_download(
  scale = 10,
  type = "rivers_lake_centerlines",
  category = "physical",
  returnclass = "sf"
)
rivers<-st_intersection(st_transform(rivers,st_crs(na)),na)
lakes <- ne_download(
  scale = 10,
  type = "lakes",
  category = "physical",
  returnclass = "sf"
)
lakes<-st_intersection(st_transform(lakes,st_crs(na)),na)
par(mar=c(0,0,0,0))
plot(st_geometry(na))
plot(st_geometry(rivers),col="blue",add=TRUE)
