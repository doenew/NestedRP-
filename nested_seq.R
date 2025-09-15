# Nested representative points based on Hybrid energy distance criterion
# nv=(n_1,n_2,\dots,n_K) # K layers and the number of points in each layer
# ts: training sample y_1,...y_N
# T: maximum number of iterations
# ifparallel: TRUE means to parallel update points in MM iterations

NestedRP_seq=function(nv,ts,T=200,ifparallel=TRUE)
{ 
  K=length(nv) # K layers
  neachlayer=c(nv[1],nv[2:K]-nv[1:(K-1)])
  z=rep(1:K,times=neachlayer) # It means which layer is each point added
  # Initial points
  N <- nrow(ts) # the number of points in training sample
  p <- ncol(ts) # dimension
  Pnc=NULL
  
  for(k in 1:K){
    Pnk= Seq_basic(Pnc =Pnc, nk=neachlayer[k],ts = ts,T = T,ifparallel =ifparallel )
    Pnc=rbind(Pnc,Pnk)
  }

  return(list(PnK=Pnc,layer=z))
}

library(randtoolbox)
nv=c(5,10,20)
Dseqmy=NestedRP_seq(nv = nv,ts=sobol(100000,2,scrambling = 3),T = 200,ifparallel = T)
#layer1

plot(Dseqmy$PnK[1:nv[1],],col=Dseqmy$layer[1:nv[1]],pch=Dseqmy$layer[1:nv[1]]+15,cex=2,xlim=c(0,1),ylim=c(0,1))
points(Dseqmy$PnK[(nv[1]+1):nv[2],],col=Dseqmy$layer[(nv[1]+1):nv[2]],pch=Dseqmy$layer[(nv[1]+1):nv[2]]+15,cex=2)
points(Dseqmy$PnK[(nv[2]+1):nv[3],],col=Dseqmy$layer[(nv[2]+1):nv[3]],pch=Dseqmy$layer[(nv[2]+1):nv[3]]+15,cex=2)
min(dist(Dseqmy$PnK))