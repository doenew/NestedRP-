

Seq_basic=function(Pnc=NULL,nk,ts,T=200,ifparallel=TRUE)
{
  # Initial augmented points
  N <- nrow(ts) # the number of points in training sample
  p <- ncol(ts) # dimension
  if(ifparallel==TRUE){
    if (!require("foreach")) install.packages("foreach")
    if (!require("doParallel")) install.packages("doParallel")
    cl <- makeCluster(detectCores()-1) # Adjust the number of threads according to the actual situation
    registerDoParallel(cl)
    cat(paste0("Parallel computing has been enabled, using ", getDoParWorkers(), " worker threads\n"))
  }
  trts=t(ts)  #The transposition of training samples is helpful for subsequent calculations
  Pnk <- ts[sample(1:N,nk,replace = FALSE),]
  Pnnew <- Pnk
  if(is.null(Pnc)){
    #MM algorithm iteration
    for(t in 1:T){
      if(ifparallel==TRUE){
        Pnnew_list <- foreach(i = 1:nk, .combine = rbind) %dopar% {
          xdy <- Pnk[i,]-trts      # Result are saved by row
          normxy <- sqrt(colSums(xdy*xdy))
          normxy[normxy==0] <- 1  # Handle zero norms
          xdx <- t(Pnk[i,]-t(Pnk))
          normxx <- sqrt(rowSums(xdx*xdx))
          normxx[normxx==0] <- 1  # Handle zero norms
          xdx=xdx/normxx  # Normalization
          
          # Update point
          update <- 1/mean(1/normxy)*(colMeans(ts/normxy)+colMeans(xdx))
          as.numeric(update)  # return a vector
          
        }
        Pnnew <- as.matrix(Pnnew_list)  
      }else{
        for(i in 1:nk){xdy <- Pnk[i,]-trts      # Result are saved by row
        normxy <- sqrt(colSums(xdy*xdy))
        normxy[normxy==0] <- 1  # Handle zero norms
        xdx <- t(Pnk[i,]-t(Pnk))
        normxx <- sqrt(rowSums(xdx*xdx))
        normxx[normxx==0] <- 1  # Handle zero norms
        xdx=xdx/normxx  # Normalization
        Pnnew[i,]<- 1/mean(1/normxy)*(colMeans(ts/normxy)+colMeans(xdx))
        
        }
         
      }
      Pnk <- Pnnew 
      cat(sprintf("迭代 %d/%d 完成\n", t, T))
      
    }
    
  }else{
    #MM algorithm iteration
    for(t in 1:T){
      if(ifparallel==TRUE){
        Pnnew_list <- foreach(i = 1:nk, .combine = rbind) %dopar% {
          xdy <- Pnk[i,]-trts      # Result are saved by row
          normxy <- sqrt(colSums(xdy*xdy))
          normxy[normxy==0] <- 1  # Handle zero norms
          xdx <- t(Pnk[i,]-t(rbind(Pnc,Pnk)))
          normxx <- sqrt(rowSums(xdx*xdx))
          normxx[normxx==0] <- 1  # Handle zero norms
          xdx=xdx/normxx  # Normalization
          
          # Update point
          update <- 1/mean(1/normxy)*(colMeans(ts/normxy)+colMeans(xdx))
          as.numeric(update)  
        }
        Pnnew <- as.matrix(Pnnew_list)  
      }else{
        for(i in 1:nk){xdy <- Pnk[i,]-trts      # Result are saved by row
        normxy <- sqrt(colSums(xdy*xdy))
        normxy[normxy==0] <- 1  # Handle zero norms
        xdx <- t(Pnk[i,]-t(rbind(Pnc,Pnk)))
        normxx <- sqrt(rowSums(xdx*xdx))
        normxx[normxx==0] <- 1  # Handle zero norms
        xdx=xdx/normxx  # Normalization
        Pnnew[i,]<- 1/mean(1/normxy)*(colMeans(ts/normxy)+colMeans(xdx))
        
        }
        
      }
      Pnk <- Pnnew 
      cat(sprintf("迭代 %d/%d 完成\n", t, T))
      
    }
    
  }
  
  return(Pnk=Pnk)
  
}















