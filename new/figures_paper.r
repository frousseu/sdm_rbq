
################################################################
################################################################
### Detailed results ###########################################
################################################################

library(mgcv)

cexright=1.1
cexrightlab=1.3
cexleft=0.9

### I vs groups
group<-c("fname")
df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
df<-merge(df,unique(d[,c("species",..group)]))
ed2<-merge(ed,df[,c("species",group)])
ed2<-ed2[order(match(ed2$species,ed$species)),]
df[,group]<-factor(df[,group],levels=rev(unique(ed2[[group]])))
df<-df[order(match(df$species,ed$species)),]
#df<-df[order(df[,group]),]
png("/data/sdm_rbq/graphics/groups.png",width=5,height=12,res=300,units="in")
par(mar=c(2,13,0.5,0.5))
form<-as.formula(paste0("I~",group))
if(group=="fname"){ex<-2}else{ex<-0}
boxplot(form,data=df[which(df$reach>=0.00),],las=2,outline=FALSE,pars=list(lwd=0.01,medlwd=1.5,medcol="grey50"),xaxt="n",yaxt="n",xlab="",ylab="",horizontal=TRUE,xaxs="i",xlim=c(1+ex,nlevels(df[,group])-ex))
for(i in pretty(0:1)){
    lines(c(i,i),c(-500,500),xpd=FALSE,lty=3,col="grey80",lwd=1)
}
for(i in 1:nlevels(df[,group])){
    lines(c(0,1),c(i,i),xpd=FALSE,lty=3,col="grey80",lwd=1)
}
points(df$I,jitter(as.integer(df[,group]),fac=if(group=="fname"){0}else{0.75}),cex=0.6,col=adjustcolor("seagreen",0.75),pch=21,lwd=0.75,bg=adjustcolor("seagreen",0.25))
at<-1:nlevels(df[,group])
labels<-levels(df[,group])
axis(2,at=at,labels=labels,mgp=c(2,0.5,0),tcl=-0.2,cex.axis=cexleft,las=2,lwd=0,lwd.ticks=1)
axis(1,mgp=c(2,0,0),tcl=-0.2,cex.axis=cexleft,lwd=0,lwd.ticks=1)
#axis(3,mgp=c(2,0.15,0),tcl=-0.2,cex.axis=0.75,lwd=0,lwd.ticks=1)
box(col="grey90",which="plot",lwd=5)
mtext(side=1,line=1.0,text="Overlap (I)")
dev.off()



### n vs reach
png("/data/sdm_rbq/graphics/I_vs_n.png",width=6,height=5,units="in",res=300)
par(mar=c(2.5,2.75,0.5,0.5))
df$logn<-log(df$n)
plot(df$n,df$I,log="x",cex=1,col=adjustcolor("seagreen",0.75),pch=21,lwd=0.75,bg=adjustcolor("seagreen",0.25),axes=FALSE,xlab="",ylab="",xlim=c(3,49000))
#abline(lm(I~n,data=df),log="x")
m<-gam(I~s(logn),family=betar(link="logit"),data=df)
#m<-loess(I~reach,data=df)
xvals<-seq(min(df$logn,na.rm=TRUE),max(df$logn,na.rm=TRUE),length.out=100)
p<-predict(m,data.frame(logn=xvals),se=TRUE,type="response")
lines(exp(xvals),p$fit,lwd=2)
polygon(exp(c(xvals,rev(xvals),xvals[1])),c(p$fit-p$se.fit*2,rev(p$fit+p$se.fit*2),(p$fit-p$se.fit*2)[1]),col=gray(0,0.1),border=NA)
grid(nx=NA,col="grey70",lty=3,lwd=2)
for(i in c(5,10,50,100,500,100,5000,10000,50000)){
    lines(c(i,i),c(-10,10),col="grey70",lty=3,lwd=2)
}
for(i in pretty(df$I)){
    lines(c(2,10^10),c(i,i),col="grey70",lty=3,lwd=2)
}
box(which="plot",lwd=6,col=gray(0,0.1))
axis(1,mgp=c(2,0.15,0),tcl=-0.2,lwd=0,lwd.ticks=1,cex.axis=cexright)
axis(2,mgp=c(2,0.25,0),tcl=-0.2,lwd=0,lwd.ticks=1,las=2,cex.axis=cexright)
mtext(side=1,line=1.5,text="Number of observations",cex=cexrightlab)
mtext(side=2,line=1.75,text="Overlap (I)",cex=cexrightlab)
dev.off()




### I vs reach
png("/data/sdm_rbq/graphics/I_vs_reach.png",width=6,height=5,units="in",res=300)
par(mar=c(2.5,2.75,0.5,0.5))
plot(df$reach,df$I,log="",cex=1,col=adjustcolor("seagreen",0.75),pch=21,lwd=0.75,bg=adjustcolor("seagreen",0.25),axes=FALSE,xlab="",ylab="")
#abline(lm(I~reach,data=df))
m<-gam(I~s(reach,k=20),data=df)
#m<-loess(I~reach,data=df)
xvals<-seq(min(df$reach,na.rm=TRUE),max(df$reach,na.rm=TRUE),length.out=100)
p<-predict(m,data.frame(reach=xvals),se=TRUE)
lines(xvals,p$fit,lwd=2)
polygon(c(xvals,rev(xvals),xvals[1]),c(p$fit-p$se.fit*2,rev(p$fit+p$se.fit*2),(p$fit-p$se.fit*2)[1]),col=gray(0,0.1),border=NA)
grid(col="grey70",lty=3,lwd=2)
box(which="plot",lwd=6,col=gray(0,0.1))
axis(1,mgp=c(2,0.15,0),tcl=-0.2,lwd=0,lwd.ticks=1,cex.axis=cexright)
axis(2,mgp=c(2,0.25,0),tcl=-0.2,lwd=0,lwd.ticks=1,las=2,cex.axis=cexright)
mtext(side=1,line=1.5,text="Proportion of asymptote reached by convex hull area",cex=cexrightlab)
mtext(side=2,line=1.75,text="Overlap (I)",cex=cexrightlab)
dev.off()

### add convexhull size

### nicheOverlap
#cl<-makeCluster(2)
#registerDoParallel(cl)
#hullarea<-foreach(i=df$species,.packages=c("sf","FNN")) %do% {
#  print(i)  
#  getobs(i) |>
#    st_union() |>
#    st_convex_hull() |>
#    st_area() |>
#    as.numeric()  
#}
#df$hullarea<-unlist(hullarea)
#stopCluster(cl)
#no<-unlist(no)

### range size vs I
png("/data/sdm_rbq/graphics/I_vs_range.png",width=6,height=5,units="in",res=300)
par(mar=c(2.5,2.75,0.5,0.5))
df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
#df$hullarea<-log(df$hullarea)
df$hullarea<-df$hullarea/(10^7)
plot(df$hullarea,df$I,log="",cex=1,col=adjustcolor("seagreen",0.75),pch=21,lwd=0.75,bg=adjustcolor("seagreen",0.25),axes=FALSE,xlab="",ylab="")
#w<-which.max(df$hullarea)
#points(df$hullarea[w],df$I[w],pch=16,col="red",cex=3)
m<-gam(I~s(hullarea),family=betar(link="logit"),data=df[which(df$reach>=0.0),])
xvals<-seq(min(df$hullarea,na.rm=TRUE),max(df$hullarea,na.rm=TRUE),length.out=100)
p<-predict(m,data.frame(hullarea=xvals),se=TRUE,type="link")
lines(xvals,boot::inv.logit(p$fit),lwd=2)
polygon(c(xvals,rev(xvals),xvals[1]),boot::inv.logit(c(p$fit-p$se.fit*2,rev(p$fit+p$se.fit*2),(p$fit-p$se.fit*2)[1])),col=gray(0,0.1),border=NA)
grid(col="grey70",lty=3,lwd=2)
box(which="plot",lwd=6,col=gray(0,0.1))
axis(1,mgp=c(2,0.15,0),tcl=-0.2,lwd=0,lwd.ticks=1,cex.axis=cexright)
axis(2,mgp=c(2,0.25,0),tcl=-0.2,lwd=0,lwd.ticks=1,las=2,cex.axis=cexright)
mtext(side=1,line=1.75,text=expression("Range size ( convex hull area "%*%~10^7~km^2~")"),cex=cexrightlab)
mtext(side=2,line=1.75,text="Overlap (I)",cex=cexrightlab)
dev.off()
### assemble
ims<-do.call("c",lapply(list.files("/data/sdm_rbq/graphics",pattern="_vs_",full=TRUE)[c(1,3,2)],image_read))
ims<-image_append(ims,stack=TRUE)
im<-image_read("/data/sdm_rbq/graphics/groups.png")
ims<-image_scale(ims,paste0("x",image_info(im)$height))
im<-image_append(c(im,ims))
im<-image_trim(im)
im<-image_border(im,"white","20x20")
#plot(im)
image_write(im,"/data/sdm_rbq/graphics/results_detailed.png")