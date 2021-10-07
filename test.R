a<-find_anchorpoints(d,30)
r<-rawmodel(y,d,a$anchorpoints,c(200,200,0.01,100),10000,8,8,"gaussian","esponenziale")
aa<-as.data.frame(a$anchorpoints)
colnames(aa)<-c("latitude","longitude")
s<-r$solutions
s<-as.data.frame(s)
colnames(s)<-c("lambda1", "lambda2", "phi", "sigma")
g<-cbind(aa,s)
#normal plot
x11()
plot(g[,1:2],pch=20)
for(i in 1:dim(aa)[1])
  draw.ellipse(x=g[i,1],y=g[i,2],a=g[i,3],b=g[i,4],angle=g[i,5],deg=FALSE)
grid()
#ggplot2
x11()
ggplot(aa, aes(x=latitude, y=longitude)) + geom_point(colour="black") + geom_ellipse(aes(x0 = latitude, y0 = longitude, a = lambda1, b = lambda2, angle = phi), data = g) + coord_fixed()
