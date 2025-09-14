# Nested representative points based on Hybrid energy distance criterion
# nv=(n_1,n_2,\dots,n_K) # K layers and the number of points in each layer
# ts: training sample y_1,...y_N
# T: maximum number of iterations
# ifparallel: TRUE means to parallel update points in MM iterations

NestedRP=function(nv,ts,T=200,ifparallel=TRUE)
{ 
  K=length(nv) # K layers
  neachlayer=c(nv[1],nv[2:K]-nv[1:(K-1)])
  z=rep(1:K,times=neachlayer) # It means which layer is each point added
  # Initial points
  N <- nrow(ts) # the number of points in training sample
  p <- ncol(ts) # dimension
  n <- nv[K] # The number of points in the maximum point set
  nstar <- sum(nv)
  Pn <- ts[sample(1:N,n,replace = FALSE),]
  
  if(ifparallel==TRUE){
    if (!require("foreach")) install.packages("foreach")
    if (!require("doParallel")) install.packages("doParallel")
    cl <- makeCluster(detectCores()-1) # Adjust the number of threads according to the actual situation
    registerDoParallel(cl)
    cat(paste0("Parallel computing has been enabled, using ", getDoParWorkers(), " worker threads\n"))
  }
  trts=t(ts)  #The transposition of training samples is helpful for subsequent calculations
  
  #MM algorithm iteration
  Pnnew <- Pn
  for(t in 1:T){
    if(ifparallel==TRUE){
      Pnnew_list <- foreach(i = 1:n, .combine = rbind) %dopar% {
        xdy <- Pn[i,]-trts      # Result are saved by row
        normxy <- sqrt(colSums(xdy*xdy))
        normxy[normxy==0] <- 1  # Handle zero norms
        xdx <- t(Pn[i,]-t(Pn))
        normxx <- sqrt(rowSums(xdx*xdx))
        normxx[normxx==0] <- 1  # Handle zero norms
        xdx=xdx/normxx  # Normalization
        cum_row_sums <- apply(xdx, 2, cumsum)
        total_avg <- numeric(ncol(xdx))
        niv <- nv[z[i]:K]
        for (i in seq_along(niv)) {
          n_i <- niv[i]
          avg_vec <- cum_row_sums[n_i, ] / n_i
          total_avg <- total_avg + avg_vec
        }
    
        final_avg <- total_avg / length(niv)

        # Update point
        update <- 1/mean(1/normxy)*(colMeans(ts/normxy)+ final_avg)
        
      }
      Pnnew <- as.matrix(Pnnew_list)  
    } else {
      for(i in 1:n){
        xdy <- Pn[i,]-trts      # Result are saved by row
        normxy <- sqrt(colSums(xdy*xdy))
        normxy[normxy==0] <- 1  # Handle zero norms
        xdx <- t(Pn[i,]-t(Pn))
        normxx <- sqrt(rowSums(xdx*xdx))
        normxx[normxx==0] <- 1  # Handle zero norms
        xdx=xdx/normxx  # Normalization
        cum_row_sums <- apply(xdx, 2, cumsum)
        total_avg <- numeric(ncol(xdx))
        niv <- nv[z[i]:K]
        for (i in seq_along(niv)) {
          n_i <- niv[i]
          avg_vec <- cum_row_sums[n_i, ] / n_i
          total_avg <- total_avg + avg_vec
        }
        
        final_avg <- total_avg / length(niv)
        
        # Update point
        update <- 1/mean(1/normxy)*(colMeans(ts/normxy)+ final_avg)
      }
    }
    Pn <- Pnnew 
    cat(sprintf("Iteration %d/%d complete\n", t, T))
  }
  
  # Close parallel cluster
  if(ifparallel==TRUE){
    stopCluster(cl)
    registerDoSEQ()  # Restore to serial computing
  }
  
  return(list(PnK=Pn,layer=z))
}


library(randtoolbox)
nv=c(10,20,30)
 D=NestedRP(nv = nv,ts=sobol(100000,2,scrambling = 3),T = 200,ifparallel = T)
 #layer1
 
 plot(D$PnK[1:nv[1],],col=D$layer[1:nv[1]],pch=D$layer[1:nv[1]]+15,cex=2,xlim=c(0,1),ylim=c(0,1))
 points(D$PnK[(nv[1]+1):nv[2],],col=D$layer[(nv[1]+1):nv[2]],pch=D$layer[(nv[1]+1):nv[2]]+15,cex=2)
 points(D$PnK[(nv[2]+1):nv[3],],col=D$layer[(nv[2]+1):nv[3]],pch=D$layer[(nv[2]+1):nv[3]]+15,cex=2)

library(support)
 Dseq1=sp(n = nv[1],p=2,dist.samp =sobol(100000,2,scrambling = 3))$sp
Dseq2=rbind( Dseq1,sp_seq(Dseq1,nv[2]-nv[1],dist.samp =sobol(100000,2,scrambling = 3))$seq)
Dseq3=rbind( Dseq2,sp_seq(Dseq2,nv[3]-nv[2],dist.samp =sobol(100000,2,scrambling = 3))$seq)

plot(Dseq3[1:nv[1],],col=D$layer[1:nv[1]],pch=D$layer[1:nv[1]]+15,cex=2,xlim=c(0,1),ylim=c(0,1))
points(Dseq3[(nv[1]+1):nv[2],],col=D$layer[(nv[1]+1):nv[2]],pch=D$layer[(nv[1]+1):nv[2]]+15,cex=2)
points(Dseq3[(nv[2]+1):nv[3],],col=D$layer[(nv[2]+1):nv[3]],pch=D$layer[(nv[2]+1):nv[3]]+15,cex=2)






