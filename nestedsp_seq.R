NestedSP_seq=function(nv,ts)
{ 
  K=length(nv) # K layers
  neachlayer=c(nv[1],nv[2:K]-nv[1:(K-1)])
  z=rep(1:K,times=neachlayer) # It means which layer is each point added
  # Initial points
  N <- nrow(ts) # the number of points in training sample
  p <- ncol(ts) # dimension
 
  library(support)
  Pnc=sp(neachlayer[1],p,dist.samp = ts)$sp
  for(k in 2:K){
    Pnk= sp_seq(D =Pnc,nseq = neachlayer[k],dist.samp = ts )$seq
    Pnc=rbind(Pnc,Pnk)
  }

  return(list(PnK=Pnc,layer=z))
}
library(randtoolbox)
nv=c(5,10,20)
Dseqsp=NestedSP_seq(nv = nv,ts=sobol(100000,2,scrambling = 3))
#layer1

plot(Dseqsp$PnK[1:nv[1],],col=Dseqsp$layer[1:nv[1]],pch=Dseqsp$layer[1:nv[1]]+15,cex=2,xlim=c(0,1),ylim=c(0,1))
points(Dseqsp$PnK[(nv[1]+1):nv[2],],col=Dseqsp$layer[(nv[1]+1):nv[2]],pch=Dseqsp$layer[(nv[1]+1):nv[2]]+15,cex=2)
points(Dseqsp$PnK[(nv[2]+1):nv[3],],col=Dseqsp$layer[(nv[2]+1):nv[3]],pch=Dseqsp$layer[(nv[2]+1):nv[3]]+15,cex=2)
min(dist(Dseqsp$PnK))