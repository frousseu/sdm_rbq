
library(mgcv)
plot(df$reach,df$I,log="",pch=16,col=gray(0,0.5))
#abline(lm(I~reach,data=df))
m<-gam(I~s(reach,k=20),data=df)
#m<-loess(I~reach,data=df)
xvals<-seq(0,max(df$reach,na.rm=TRUE),length.out=100)
p<-predict(m,data.frame(reach=xvals),se=TRUE)
lines(xvals,p$fit,lwd=2)
polygon(c(xvals,rev(xvals),xvals[1]),c(p$fit-p$se.fit*2,rev(p$fit+p$se.fit*2),(p$fit-p$se.fit*2)[1]),col=gray(0,0.1),border=NA)

df$logn<-log(df$n)
plot(df$n,df$I,log="x",pch=16,col=gray(0,0.5))
#abline(lm(I~n,data=df),log="x")
m<-gam(I~s(logn),family=betar(link="logit"),data=df)
#m<-loess(I~reach,data=df)
xvals<-seq(min(df$logn,na.rm=TRUE),max(df$logn,na.rm=TRUE),length.out=100)
p<-predict(m,data.frame(logn=xvals),se=TRUE,type="response")
lines(exp(xvals),p$fit,lwd=2)
polygon(exp(c(xvals,rev(xvals),xvals[1])),c(p$fit-p$se.fit*2,rev(p$fit+p$se.fit*2),(p$fit-p$se.fit*2)[1]),col=gray(0,0.1),border=NA)



df<-read.csv("/data/sdm_rbq/graphics/mapSpeciesres.csv")
df<-merge(df,unique(d[,c("species","fname")]))
ed2<-merge(ed,df[,c("species","fname")])
ed2<-ed2[order(match(ed2$species,ed$species)),]
df$fname<-factor(df$fname,levels=rev(unique(ed2$fname)))


png("/data/sdm_rbq/graphics/groups.png",width=5,height=12,res=300,units="in")
par(mar=c(1,11,1,0.5))
boxplot(I~fname,data=df,las=2,outline=FALSE,pars=list(lwd=0.01,medlwd=1.5,medcol="grey50"),xaxt="n",yaxt="n",xlab="",ylab="",horizontal=TRUE,xaxs="i",xlim=c(1+2,nlevels(df$fname)-2))
for(i in pretty(0:1)){
    lines(c(i,i),c(-500,500),xpd=FALSE,lty=3,col="grey80",lwd=1)
}
for(i in 1:nlevels(df$fname)){
    lines(c(0,1),c(i,i),xpd=FALSE,lty=3,col="grey80",lwd=1)
}
points(df$I,as.integer(df$fname),cex=0.6,col=adjustcolor("seagreen",0.75),pch=21,lwd=0.75,bg=adjustcolor("seagreen",0.25),bty="n")
at<-1:nlevels(df$fname)
labels<-levels(df$fname)
axis(2,at=at,labels=labels,mgp=c(2,0.5,0),tcl=-0.2,cex.axis=0.75,las=2,lwd=0,lwd.ticks=1)
axis(1,mgp=c(2,0,0),tcl=-0.2,cex.axis=0.75,lwd=0,lwd.ticks=1)
axis(3,mgp=c(2,0.15,0),tcl=-0.2,cex.axis=0.75,lwd=0,lwd.ticks=1)
box(col="grey90",which="plot",lwd=5)
dev.off()


