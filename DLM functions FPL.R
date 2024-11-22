# Functions to run the EM algorithm, optimize discount factor and run DLM

update.mu0 <- function(mu, Gt, start.time){
  
  for(i in 1:(start.time-1)){
    mu <- Gt %*% mu
  }
  
  return(mu)
  
}


## Function to estimate the mu0 (initial parameter vector) ----
get.mu0 <- function(Data, stratify.by=NA, time.var, expected.start.time, 
                    relevant.names, simple.linear = FALSE){
  
  D.full <- Data
  mu0.list <- list()
  
  # If stratify.by is NA, we can handle it by adding a dummy column
  if(is.na(stratify.by)){
    D.full$Dummy <- 1
    stratify.by <- 'Dummy'
  }
  for(stratification.group in sort(unique(D.full[,stratify.by]))){
    
    # Prepare for making mu0
    Data <- subset(D.full, D.full[,stratify.by] == stratification.group)
    mu0 <- c()
    mu0.names <- c()
    
    for(name in relevant.names){ 
      
      # Elements for mu0
      D.B <- subset(Data, !is.na(Data[,name]) & !is.na(Data[,time.var]))
      x <- D.B[,time.var]
      y <- D.B[,name]
      Spline = smooth.spline(x = x, y = y, df = 4)
      
      # - plot to see what's going on
      # plot(y~x, xlab=time.var, ylab=name, main = paste(stratify.by, '=', stratification.group))
      # lines(Spline, col='red', lwd=3)
      # plot(Spline, type='l', lwd=2, ylab=name, xlab=time.var, 
      #      main = paste(stratify.by, '=', stratification.group))
      
      # Finalize m0
      pred0 <- predict(object = Spline, x = (expected.start.time-1))$y
      
      if(simple.linear == FALSE){
        mu0 <- c(mu0, c(pred0, 1))
      }else{
        pred1 <- predict(object = Spline, x = (expected.start.time))$y
        Diff <- pred1-pred0
        mu0 <- c(mu0, c(pred0, Diff))
      }
      
      mu0.names <- c(mu0.names,
                     name,
                     paste('d.', name, sep=''))
      
      mu0 <- matrix(mu0)
      row.names(mu0) <- mu0.names
    }
    
    # Save the estimated parameter vector
    mu0.name <- paste('mu0_', stratification.group, sep='')
    mu0.list[[length(mu0.list)+1]] <- mu0
    names(mu0.list)[length(mu0.list)] <- mu0.name
  }
  
  return(mu0.list)
}


## Function to estimate the C0 (initial prior variance) ----
get.C0 <- function(Data, stratify.by,
                   time.var,
                   expected.start.time,
                   relevant.names){
  
  D.full <- Data
  C0.list <- list()
  
  # If stratify.by is NA, we can handle it by adding a dummy column
  if(is.na(stratify.by)){
    D.full$Dummy <- 1
    stratify.by <- 'Dummy'
  }
  for(stratification.group in sort(unique(D.full[,stratify.by]))){
    
    # Prepare for making C0
    D <- subset(D.full, D.full[,stratify.by] == stratification.group)
    Data.A <- D
    C0.names <- c()
    
    for(name in relevant.names){
      
      # Elements for C0
      # - diff for the prior variance on the initial rate of change in mean weight
      Diff <- diff(D[,name])
      Data.A[,paste('d.', name, sep='')] <- c(NA,Diff)
      
      C0.names <- c(C0.names,
                    name,
                    paste('d.', name, sep=''))
    }
    
    # Make C0
    # - Re-order the colums
    # Data.A <- Data.A[,c(metadata.names, mu0.names)]
    Data.A <- Data.A[,c(time.var, C0.names)]
    
    # - Only consider the start of the lactation for C0
    # - omit the actual first values per time series (expected.start.time),
    #   as these will have extreme outliers in the Diff values, 
    #   stemming from the differences between the first vlaue in a new series 
    #   and the last value en the previous series
    Data.A <- subset(Data.A, Data.A[,time.var] %in% (expected.start.time+1):(expected.start.time+2))
    Data.A <- na.omit(Data.A)
    C0 <- cov(Data.A[,C0.names])
    rm(Data.A)
    
    # Save the estimated variance component
    C0.name <- paste('C0_', stratification.group, sep='')
    C0.list[[length(C0.list)+1]] <- C0
    names(C0.list)[length(C0.list)] <- C0.name
  }
  
  return(C0.list)
}


## Function for getting spline functions to be used in Gt matrix ----
get.spline <- function(Data,
                       stratify.by=NA,
                       time.var,
                       relevant.names,
                       plot.it = FALSE){

  # If stratify.by is na, we need to add a Dummy variabe
  if(is.na(stratify.by)){
    Data$Dummy <- 1
    stratify.by <- 'Dummy'
  }
  
  D.full <- Data
  Spline.list <- list()
  
  for(name in relevant.names){ 

    spline.list.name <- list()
    counter <- 0
    for.limits <- c()
    for(stratification.group in sort(unique(D.full[,stratify.by]))){
      
      counter <- counter + 1
      
      print(paste(name, '|', stratification.group))

      Data <- subset(D.full, D.full[,stratify.by] == stratification.group)
      
      # Create the spline
      D.B <- subset(Data, !is.na(Data[,name]))
   
      if(dim(D.B)[1] > 0){
        x <- D.B[,time.var]
        y <- D.B[,name]
        
        agg <- aggregate(x = D.B[,name], by=list(D.B[,time.var]), FUN=mean)
     
        # Make a spline function with 4 degrees of freedom
        Spline = smooth.spline(x = x, y = y, df = 4)
        for.limits <- c(for.limits, range(Spline$y))
        
        for.limits <- c(for.limits, range(c(Spline$y)))
        
        # Save the spline function
        spline.name <- paste(name, '_Spline_', stratification.group, sep='')
        Spline.list[[length(Spline.list)+1]] <- Spline
        names(Spline.list)[length(Spline.list)] <- spline.name
        
        # A secondary list, just for plotting
        if(plot.it == TRUE){
          spline.list.name[[length(spline.list.name)+1]] <- Spline
        }
      }else{
        print(paste('ERROR! No non-missing data for', name, '&', stratification.group))
      }
      
    }
    
    # Plot it, if relevant
    if(plot.it == TRUE){
      for(i in 1:length(spline.list.name)){
        Spline <- spline.list.name[[i]]
        if(i == 1){
          ylim=c(range(c(Spline$y, agg$x)))
          plot(Spline$y~Spline$x, xlab=time.var, ylab=name, ylim=ylim, type='l', lty=i, main=relevant.names)
          points(agg)
        }else{
          lines(Spline$y~Spline$x, lty=i)
        }
      }
    }
    
  }
  return(Spline.list)
}


## Make Gt matrix
get.Gt <- function(Data.A=Data_,
                   i.A=i,
                   time.var.A=time.var,
                   stratify.by.A=stratify.by,
                   Spline.list.A=Spline.list,
                   relevant.names.A=relevant.names){
  
  Data_ <- Data.A
  i <- i.A
  time.var <- time.var.A
  stratify.by <- stratify.by.A
  Spline.list <- Spline.list.A
  relevant.names <- relevant.names.A
  
  na.spline.i <- which(is.na(spline.list))
  
  # If stratify.by is na, we need to add a Dummy variabe
  if(is.na(stratify.by)){
    Data_$Dummy <- 1
    stratify.by <- 'Dummy'
  }
  
  # Make the diagonal of the Gt matrix
  Gt <- diag(1,length(relevant.names)*2) 

  # Get the relevant trend component for each variable for this observation
  for(name in relevant.names){
    
    # We expect all parameters to follow a 2x2 structure
    k <- which(relevant.names == name)
    i.A <- 2*k - 1
    j <- 2*k
    
    # - get the relevant spline function
    stratification.group <- Data_[i,stratify.by]
    spline.name <- paste(name, '_Spline_', stratification.group, sep='')
    spline.name <<- spline.name
    Spline <- Spline.list[[which(names(Spline.list)==spline.name)]]
    
    # - get the relevant spline function for the upper-right element
    if(identical(Spline.list, NA)){
      Trend <- 1
    }else{

      # - check if the spline function for the current name is NA - then we want no trend!
      if(identical(Spline, NA)){
        Trend <- 0
      }else{
        time <- Data_[i,time.var]
        time_1 <- Data_[i-1,time.var]
        pred <- predict(object = Spline, x = time)$y
        pred_1 <- predict(object = Spline, x = time_1)$y
        if(length(pred_1) == 0){
          pred_1 <- 0
        }
        Trend <- pred - pred_1
        }
      
    }
   
    Gt[i.A,j] <- Trend
    
  }
  
  return(Gt)
}


## Make Ft matrix
get.Ft <- function(relevant.names.A = relevant.names){
  
  relevant.names <- relevant.names.A
  
  N <- (length(relevant.names)*2) * length(relevant.names)
  Ft <- matrix(rep(0,N), nrow = length(relevant.names))
  for(name in relevant.names){
    k <- which(relevant.names == name)
    i1 <- k
    j1 <- 2*k-1
    Ft[i1,j1] <- 1
  }
  
  Ft <- t(Ft)
  
  return(Ft)
}



# Moving average, moving SData, you name it
moving.function <- function(x, n, FUN){
  
  start <- floor(n/2)+1
  N <- length(x) - floor(n/2)
  out <- c()
  if(N > start){
    for(i in start:N){
      obs <- x[(i-floor(n/2)):(i+floor(n/2))]
      res <- FUN(na.omit(obs))
      out <- c(out,res)
    }
    NAs <- rep(NA, floor(n/2))
    out <- c(NAs, out, NAs)
  }else{
    out <- rep(NA, length(x))
  }
  
  return(out)
}


## Function to estimate V (observational variance) ----
get.V <- function(Data, identifyer=NA, stratify.by=NA, time.var, relevant.names){
  

  D.full <- Data
  V.list <- list()
  
  # If identifyer is NA, we can handle it by adding a dummy column
  if(is.na(identifyer)){
    D.full$Dummy_ID <- 1
    identifyer <- 'Dummy_ID'
  }
  # If stratify.by is NA, we can handle it by adding a dummy column
  if(is.na(stratify.by)){
    D.full$Dummy <- 1
    stratify.by <- 'Dummy'
  }
  for(stratification.group in sort(unique(D.full[,stratify.by]))){
    
    # Prepare for making V
    D <- subset(D.full, D.full[,stratify.by] == stratification.group)
    for.V <- cbind()
    
    for(name in relevant.names){
      
      # Normalize the name's data - only relevant for multivariate data
      # if(length(relevant.names) > 1){
      #   D[,name] <- D[,name]/Normalization.factors[,name]
      # }
      
      # Elements for the observational matrix, V
      resiData.All <- c()
      for(ID in unique(D[,identifyer])){
        
        ID.set <- subset(Data, Data[,identifyer] == ID)
        
        #use a two-sided moving average to estimate V
        y <- ID.set[,name]
        x <- ID.set[,time.var]
        MA <- moving.function(x = y, n = 5, FUN = mean)
        
        # plot(y~x, xlab=time.var, ylab=name)
        # lines(MA, col='red')
        
        #find the residuals between the observed and filtered/estimated values
        resid.ID <- ID.set[,name] - MA
        
        #combine the residuals for all individual pen-batches to estimate an overall V
        resiData.All <- c(resiData.All, resid.ID)
      }
      
      # Plot the histogram to make sure the residuals are normally distributed around 0
      if(length(resiData.All) > 100){
        hist(resiData.All, round(length(resiData.All)/10),
             main = paste('residuals', name))
      }else{
        hist(resiData.All)
      }
      abline(v=0, col='red')
      
      # Add the residuals to a data frame, which we will use for making the final V matrix
      a <- paste('resid.', name, sep='')
      for.V <- cbind(for.V,
                     resiData.All)
      colnames(for.V)[ncol(for.V)] <- a
      
    }
    
    # Estimate the observational variance matrix, V
    colnames(for.V) <- relevant.names
    V <- cov(na.omit(for.V))
    
    # Save the estimated variance component
    V.name <- paste('V_', stratification.group, sep='')
    V.list[[length(V.list)+1]] <- V
    names(V.list)[length(V.list)] <- V.name
  }
  
  return(V.list)
}


# Get a martrix with 1 in the cells where the observation 
# contributes to the observation variance-covariance matrix.
# Other cells are 0 (When it's NA). Used by the EM-algorithm
getVSumElement = function(Data, n) {
  row = Data[n, ]
  rem = c()
  V = matrix(1, nrow=length(relevant.names), ncol=length(relevant.names))
  for(cont.var in relevant.names){
    if (is.na(row[cont.var])) { #cont.var=relevant.names[1]
      cont.var.index <- which(relevant.names == cont.var)
      V[cont.var.index, ] = 0
      V[, cont.var.index] = 0
    }
  }
  
  return(V)
}


## A function for running the DLM
runDLM <- function(Data,
                   mu0,
                   C0,
                   V,
                   W=NA,
                   adjust.W=FALSE,
                   delta=0.95,
                   relevant.names,
                   Spline.list=NA,
                   time.var,
                   stratify.by=NA){
  
  n <- nrow(Data)
  
  Yt.list <- list()
  at.list <- list()		 # Define the lists
  Rt.list <- list()
  ft.list <- list()
  Qt.list <- list()
  At.list <- list()
  et.list <- list()
  ut.list <- list()
  mt.list <- list()
  Ct.list <- list()
  Ft.list <- list()
  VSE.list <- list()
  Gt.list <- list()
  

  mt <- mu0				# Prior Distribution
  Ct <- C0
  # Make sure Ct is symmetrical
  Ct <- (Ct + t(Ct))/2
  Ct <<- Ct
  
  #When we have an multivariate DLM with trend, both mu0 and C0 are not numbers, 
  #but matrices. Thereby, remember to apply %*% in all cases. 
  
  for(i in (1:n)){
    
    # Define relevant variables as global variables
    # - it's not pretty, but it makes it easier to use custom get.Gt and get.Ft functions
    Data_ <<- Data
    i <<- i
    time.var <<- time.var
    stratify.by <<- stratify.by
    Spline.list <<- Spline.list
    relevant.names <<- relevant.names
    
    # Get V-sum-element, to be used in the EM-algorithm
    VSE <- getVSumElement(Data, i)
    
    # Get the observation vector
    # Yt <- getYt(Data, i, relevant.names)
    Yt <- t(as.matrix(Data[i,relevant.names]))
    colnames(Yt) <- NULL
    
    # Get the observational variance (Vt), including only values related to observed variables
    # Vt <- getVt(Data, i, V, relevant.names)
    Vt <- V
    Vt <<- Vt
    
    # Define Wt
    if(identical(W, NA)){
      Wt <- ((1-delta)/delta) * Ct
      Wt <- (Wt + t(Wt))/2
    }else{
      Wt <- W
    }
    Wt <<- Wt
    
    # Make the DLM more adaptive in the beginning
    if(adjust.W == TRUE){
      
      if(i == 1){
           Wt <- Wt * 20000
      }
   

      # ADDITION SPECIFIC FOR THIS FPL MODEL:
      # Make the DLM more adaptive when type is changed
      if(Data[i, ]$Switched == 1 & !is.na(Data[i,]$Switched)){
        Wt <- Wt * 20000
      }
      
    }
  
   
    
    # Get Gt - independent of the current Yt
    Gt <- get.Gt()
    # Gt <- Gt[-spline.na.i,-spline.na.i]
    Gt <<- Gt
    
    # Get Ft given the current Yt
    # - only rows related to observed variables are included
    Ft <- get.Ft()
    # Ft <- Ft[-spline.na.i,] 
    Ft <<- Ft
    
    # Run the Kalman filter - only if we observe at least one variable!
    mt <<- mt
    at <- Gt %*%  mt		                   # Prior mean
    Rt <- Gt %*%  Ct %*% t(Gt) + Wt        #! I have changes it to G'     # Prior Variance
    Rt <<- Rt
    
    ft <- t(Ft) %*% at		      	         # One-step Forecast mean
    Qt <- t(Ft) %*% Rt %*%  Ft + Vt        # One-step Forecast variance
    
    At <- Rt %*% Ft %*% solve(Qt)          # Adaptative Coef. matrix
    et <- Yt  - ft	                       # one-step forecast error
    ut <- et / sqrt(diag(Qt))              #Standardized forecast error
    
    # - handle the missing values
    et.A <- et
    et.A[which(is.na(et.A))] <- 0
    
    # - update the parameter vector and variance matrix
    mt <- at + At %*% et.A                # Filtered mean
    Ct <- Rt - At  %*% Qt %*% t(At)	      # Filtered variance
    
    # - if there are missing values, undo the updating we just did fto Ct for the relevant variables!
    if(length(which(is.na(et)) > 0)){
      
      na.names <- rownames(et)[which(is.na(et))]
      Ct_na.i <- c()
      for(name in na.names){
        Ct_na.i <-c(Ct_na.i, grep(pattern = name, x = rownames(Ct)) )
      }
      
      Ct[Ct_na.i, Ct_na.i] <- Rt[Ct_na.i, Ct_na.i]
    }
    
    
    # Make sure Ct is symmetrical
    Ct <- (Ct + t(Ct))/2
    Ct <<- Ct
    

    # Save the values in lists
    Yt.list[[i]] <- Yt
    at.list[[i]] <- at
    Rt.list[[i]] <- Rt
    ft.list[[i]] <- ft
    Qt.list[[i]] <- Qt
    At.list[[i]] <- At
    et.list[[i]] <- et
    ut.list[[i]] <- ut
    mt.list[[i]] <- mt
    Ct.list[[i]] <- Ct
    Ft.list[[i]] <- t(Ft)
    VSE.list[[i]] <- VSE
    Gt.list[[i]] <- Gt
    
  }
  
  # rm(Data)
  # rm(i)
  # rm(time.var)
  # rm(stratify.by)
  # rm(Spline.list)
  # rm(relevant.names)
  
  return(list(
    Yt=Yt.list,
    at=at.list,
    Rt=Rt.list,
    ft=ft.list,
    Qt=Qt.list,
    At=At.list,
    et=et.list,
    ut=ut.list,
    mt=mt.list,
    Ct=Ct.list,
    F=Ft.list,
    vse=VSE.list,
    Gt.list=Gt.list
  ))
}



# The Kalman smoother applied to a D in a pen.
# The input parameter is a result returned from the filter.
runSmoother <- function(res) {
  
  n = length(res$mt)
  p = length(res$mt[[1]])
  mts <- array(NA,dim=c(p,1,n));
  Cts <- array(NA,dim=c(p,p,n));
  
  # Put last value equal to filtered
  mts[,,n] <- res$mt[[n]]
  Cts[,,n] <- res$Ct[[n]]
  
  # These are useful
  Bt <- array(NA,dim=c(p,p,n))
  Lt <- array(NA,dim=c(p,p,n));  
  
  # Iterate backwards over days
  for(i in ((n-1):1))   {
    
    # Get Gt
    Gt <- res$Gt.list[[i+1]]
    
    #print(i)
    res$R[[i+1]] <- as.matrix(res$R[[i+1]])
    
    Bt[,,i] <- as.matrix( res$Ct[[i]] %*% t(Gt) %*% solve(res$R[[i+1]]) )        #MCMC - page 570 on the book (West and Harrison, 1997)
    mts[,,i] <- res$mt[[i]] + Bt[,,i] %*% (mts[,,i+1] - res$a[[i+1]])
    Cts[,,i] <- as.matrix( res$C[[i]] + Bt[,,i] %*% (Cts[,,i+1] - res$R[[i+1]]) %*% t(Bt[,,i]) )
  }
  
  # Now when we are at it: Find L and store it for the EM algorithm
  for(i in ((n):2))  {
    Lt[,,i] <- Cts[,,i] + Gt%*%Cts[,,i-1]%*%t(Gt) - Cts[,,i]%*%t(Bt[,,i-1]) - Bt[,,i-1]%*%Cts[,,i]
  }
  
  return(list(mts=mts,
              Cts=Cts,
              Lt=Lt,
              Data=res$D));
}


## Run the EM algorithm over the indicated number of steps
runEM = function(Des, mu0, C0, V0, W0, steps = 1, silent = TRUE, 
                 DLM.version = runDLM, relevant.names, Spline.list, 
                 time.var, stratify.by) {
  Vs = list()
  Ws = list()
  Vn = V0
  Wn = W0
  Cn <- diag(0, length(mu0))
  mu0n <- rep(0, length(mu0))
  

  #Choose wich version of DLM to use
  runDLM <- DLM.version
  
  # Iterate over steps
  for (s in 1:steps) {
    print(paste("Step", s))
    
    # Set sums and counts to zero
    # Count structure for observation variance
    sumV = matrix(0, nrow=length(relevant.names), ncol=length(relevant.names))
    # Sum for observation variance
    sumObs = matrix(0, nrow=length(relevant.names), ncol=length(relevant.names))
    # Count for system variance
    sumW = 0
    # Sum for observation variance
    sumSys = matrix(0, nrow=length(mu0), ncol=length(mu0))
    # Iterate over Des
    for (b in 1:length(Des)) {
      
      # Expectation step - filter and smooth
      if(silent == FALSE){
        progress <- b/length(Des)*100
        print(paste(progress, '%'))
      }
      
     
      # Get mu0 and C0, given the first month of our data set
      Data = Des[[b]]

      # Get mu0 and C0 - NOTICE: THIS DOES NOT DO ANYTHING, SINCE MU0 AND C0 ARE DEFINED IN THE DLM
      # mu0 <- get.mu0(Data = training.set, stratify.by=NA, 
      #            time.var = time.var.names.All, expected.start.time = Data$HatchMonth[1], 
      #            relevant.names, simple.linear = FALSE)$mu0_1
      # 
      # C0 <- get.C0(Data = training.set, stratify.by=NA, 
      #          time.var = time.var.names.All, expected.start.time = Data$HatchMonth[1], 
      #          relevant.names)$C0_1
      # 
      # mu0 <- as.matrix(mu0[-spline.na.i,])
      # C0 <- C0[-spline.na.i,-spline.na.i]
      # 
    
      a <- get.mu0_and_C0(Data = learning.set,
                    relevant.names,
                    time.var.names.All,
                    trend.all, 
                    N.w.all,
                    w.all)

      Gt <- get.Gt()
      mu0 <- a$mu0
      C0 <- a$C0
      # Update mu0 to correspond with the start time of the data
      # mu0 <- update.mu0(mu0, Gt, start.time=1)
      
      # Run the DLM
      res <- runDLM(Data = Des[[b]], mu0 = mu0, C0 = C0, V = Vn, W = Wn, adjust.W=FALSE, delta=NA, relevant.names = relevant.names, Spline.list=Spline.list, time.var=time.var, stratify.by)
      
      smot = runSmoother(res)
      
      # Get the smoothened C0 from this one
      Cn <- Cn + smot$Cts[,,1]
      mu0n <-  mu0n + smot$mts[,,1]
      
      # Set contributions to sums and counts to zero for this D
      bSumV = matrix(0, nrow=length(relevant.names), ncol=length(relevant.names))
      sumCountV = matrix(0, nrow=length(relevant.names), ncol=length(relevant.names))
      bSumW = matrix(0, nrow=length(mu0), ncol=length(mu0))
      sumCountW = 0
      
      # Sum contributions
      n = length(res$mt)
      
      # Iterate over time within D 
      for (t in 1:n) {
        
        # Get the Gt 
        Gt <- res$Gt.list[[t]]
        
        # Only if observations at all
        if (length(res$Yt[[t]]) > 0) {
          
          # Observation variance
          
          # Find the contribution to the sum even though it does not have the correct dimension
          Vcont = res$F[[t]]%*%smot$Cts[,,t]%*%t(res$F[[t]]) + 
            (res$Y[[t]] - res$F[[t]]%*%smot$mts[,,t])%*%t((res$Y[[t]] - res$F[[t]]%*%smot$mts[,,t]))       #can't find this in the book 
          # Get the pointer matrix
          vse = res$vse[[t]]
          # Create a full 3 x 3 matrix
          Vfull = matrix(0, nrow = length(relevant.names), ncol = length(relevant.names))
          # Find out which cells to enter
          ind = c()
          for (i in 1:length(relevant.names)) {
            if (vse[i,i] > 0) {
              ind = c(ind, i)
            }
          }
          
          # Enter the contribution into the right cells of the 3 x 3 matrix
          Vfull[ind,ind] = Vcont[ind,ind]
          # Add the resulting 3 x 3 matrix
          bSumV = bSumV + Vfull
          # Adjust the counts matrix
          sumCountV = sumCountV + vse          
          
          
          # System variance
          if (t > 1) {
            # Find the contribution - the dimension is always correct
            bSumW = bSumW + smot$Lt[,,t] +                                                            #can't find this in the book 
              (smot$mts[,,t] - Gt%*%smot$mts[,,t-1])%*%t(smot$mts[,,t] - Gt%*%smot$mts[,,t-1])       
            sumCountW = sumCountW + 1
          }
          
        }
      }
      
      # Check for negative variances
      ignore = FALSE
      for (j in 1:length(relevant.names)) {
        if (bSumV[j, j] < 0) {
          # Adjust to 0
          bSumV[j, j] = 0
          bSumV[j, ] = 0
          bSumV[, j] = 0
          # Print a comment
          if (! silent) print(paste("Negative contribution to observation variance", j, "for D", b))
        }
      }
      
      # Check for negative variances
      ignore = FALSE
      for (j in 1:length(mu0)) {
        if (bSumW[j, j] < 0) {
          # Adjust to 0
          bSumW[j, j] = 0
          bSumW[, j] = 0
          bSumW[j, ] = 0
          # Print a message
          if (! silent) print(paste("Negative contribution to system variance", j, "for D", b))
        }
      }
      
      # This will never happen (used for debugging)
      if (ignore) {
        print(paste("Contribution to observation variance ignored from D", b))
        for (i in 1:length(relevant.names)) {
          for (j in 1:length(relevant.names)) {
            bSumV[i, j] = 0
          }
        }
        sumCountV = matrix(0, nrow=length(relevant.names), ncol=length(relevant.names))
      }
      # Add the contribution from the D to the total
      sumObs = sumObs + bSumV
      sumV = sumV + sumCountV
      
      
      # This will never happen (used for debugging)
      if (ignore) {
        print(paste("Contribution to system variance ignored from D", b))
        for (i in 1:length(mu0)) {
          for (j in 1:length(mu0)) {
            bSumW[i, j] = 0
          }
        }
        sumCountW = 0
      }
      # Add the contribution from the D to the total
      sumSys = sumSys + bSumW
      sumW = sumW + sumCountW
      
    }
    # Normalize by counts
    Vn = sumObs/sumV
    Wn = sumSys/sumW
    Cn <- Cn/n
    mu0n <-  mu0n/n
    # Make sure they are symmetric
    Vn = (Vn + t(Vn))/2
    Wn = (Wn + t(Wn))/2
    Cn <- (Cn + t(Cn))/2
    Vn <<- Vn
    Wn <<- Wn
    # Mke sure they have the right names
    colnames(Vn) <- relevant.names
    rownames(Vn) <- relevant.names
    
    colnames(Wn) <- rownames(mu0)
    rownames(Wn) <- rownames(mu0)
    
    colnames(Cn) <- rownames(mu0)
    rownames(Cn) <- rownames(mu0)
    
    mu0n <- as.matrix(mu0n)
    rownames(mu0n) <- rownames(mu0)
    
    # Save them in lists
    Vs[[s]] = Vn
    Ws[[s]] = Wn
  }
  return(list(V=Vs, W=Ws, Cn=Cn, mu0n=mu0n))
}


## Run the EM algorithm with early stopping
runEM_earlyStopping <- function(Data, stratify.by, Spline.list, identifyer, 
                                V0=NA, W0=NA, C0.list, mu0.list, 
                                no.better.limit=1, time.var=NA, relevant.names, 
                                round.by=4){
  
  ### Run the EM algorithm

  # Prepare output lists
  V.list.out <- list()
  W.list.out <- list()
  mu0.list.out <- list()
  C0.list.out <- list()
  
  # If identifyer is NA, we can handle it by adding a dummy column
  if(is.na(identifyer)){
    Data$Dummy_ID <- 1
    identifyer <- 'Dummy_ID'
  }
  # Iterate over the stratification groups
  if(is.na(stratify.by)){
    Data$Dummy <- 1
    stratify.by <- 'Dummy'
  }
  D.full <- Data
  for(stratification.group in unique(D.full[,stratify.by])){
    
    D <- subset(D.full, D.full[,stratify.by] == stratification.group)
    divide.start <- Sys.time()
    
    # Randomly select a training set and a test set to determine when to stop the EM algorithm
    set.seed(111)
    print('Dividing the data into a training and test set for EM ...')
    # test.IDs <- as.data.frame(sample(x = unique(D[,identifyer]), size = ceiling(0.3*length(as.numeric(unique(D[,identifyer])))), replace = FALSE))
    test.IDs <- sample(x = unique(D[,identifyer]), 
                       size = ceiling(0.3*length(as.numeric(unique(D[,identifyer])))), 
                       replace = FALSE)
    training.IDs <- unique(D[,identifyer])[which(unique(D[,identifyer]) %in% test.IDs == FALSE)]
    if(length(training.IDs)==0){
      training.IDs <- test.IDs
    }
    
    test.set <- D[which(D[,identifyer] %in% test.IDs),]
    training.set <<- D[which(D[,identifyer] %in% training.IDs),]
    # test.set <- subset(Data, D[,identifyer] %in% test.IDs)
    
    # - we need all the training IDs to be elements of a list (that is how the EM algorithm takes them)
    start.makeList <- Sys.time()
    Des <- list()
    for(ID in training.IDs){
      progress <- which(training.IDs == ID)/length(training.IDs)*100
      # print(paste('Making list of training IDs |', progress, '%'))
      ID.set <- subset(Data, D[,identifyer] == ID)
      Des[[length(Des)+1]] <- ID.set
    }
    # print(Sys.time()-start.makeList)
    
    train.set <- Des
    
    print(paste('--Done! It took', round(difftime(Sys.time(), divide.start, units = 'min'),1), 'minutes'))
    
    # Run the EM algorithm until it no longer improves the performance of the DLM
    print('Running the EM algorithm ...')
    

    
    # Fist we need to see if we have stratification-specific versions of C0
    Grep <- grep(pattern = stratification.group, x = names(C0.list))
    if( length(Grep) == 1 ){
      C0 <- C0.list[[Grep]]
    }else{
      # if not, we just use the first element of C0
      C0 <- C0.list[[1]]
    }
    C <- C0
    
    # Now we need to see if we have stratification-specific versions of mu0
    Grep <- grep(pattern = stratification.group, x = names(mu0.list))
    if( length(Grep) == 1 ){
      mu0 <- mu0.list[[Grep]]
    }else{
      # if not, we just use the first element of C0
      mu0 <- mu0.list[[1]]
    }
    mu <- mu0
    
    # V0 have the option of being NA - we need to handle that, so we arbitrarily use 1/10 of the relevant diagonal values from C0
    if(identical(V0, NA)){
      names.i <- which(colnames(C) %in% relevant.names)
      V0 <- C[names.i,names.i]/10
    }
    V = V0
    
    # W0 have the option of being NA - we need to handle that, so we arbitrarily use 1/10 of the values from C0
    if(identical(W0, NA)){
      names.i <- which(colnames(C) %in% relevant.names)
      W0 <- C/10
    }
    W = W0
   
    # Set initial best values and run the EM algorithm
    no.better <- 0
    RMSE.best <- 9999999999
    total.steps <- 0
    V.best <- V
    W.best <- W
    C0.best <- C
    mu0.best <- mu
    while(no.better < no.better.limit){
      
      # Run it on the training set
      start.time <- Sys.time()
      a <- runEM(Des = train.set, mu0 = mu, C0 = C, V0 = V, W0 = W, 
                 steps = 1, silent = TRUE, DLM.version = runDLM, 
                 relevant.names, Spline.list, time.var, stratify.by)
                 # Des,             mu0,      C0,     V0,     W0,     steps = 1, silent = TRUE, DLM.version = runDLM, relevant.names,                  Spline.list, time.var, stratify.by
      # print(Sys.time()-start.time)
      V <- a$V[[1]]
      W <- a$W[[1]]
      # C <- a$Cn
      # mu <- a$mu0n
      
      # Try running the DLM with these variance parameters on all PB's in the test set
      et.all <- cbind()
      Yt.all <- cbind()
      # et.all <- c()
      for(ID in unique(test.set[,identifyer])){
        test.set.i <- subset(test.set, test.set[,identifyer] == ID)
        res <- runDLM(Data = test.set.i, mu0 = mu,  C0 = C,  V=V,    W=W,    
                      adjust.W=FALSE, delta=NA, relevant.names = relevant.names, 
                      Spline.list=Spline.list, time.var=time.var, stratify.by)
       
        # - get the raw forecast errors as a data frame
        et <- t(as.data.frame(res$et))
        row.names(et) <- NULL
        et.all <- rbind(et.all, et)
        
        # - get the observations as a data frame
        Yt <- t(as.data.frame(res$Yt))
        row.names(Yt) <- NULL
        Yt.all <- rbind(Yt.all, Yt)
        
        # # -- as a vector
        # et.all <- c(et.all, unlist(res$et))
      }
      et.all <- as.data.frame(et.all)
      row.names(et.all) <- NULL
      
      # # # Normalize them for a fair comparison
      # # par(mfrow=c(3,2))
      # if(length(relevant.names) > 1){
      #   for(i in 1:ncol(et.all)){
      #     # hist(et.all[,i],100, main=relevant.names[i])
      #     # et.all[,i] <- et.all[,i]/Means[i]
      #     et.all[,i] <- et.all[,i]/(Yt.all[,i]+1)
      #     # hist(et.all[,i],100, main=relevant.names[i])
      #   }
      # }
      
      
      # Calculate the RMSE and check if it is still improving
      RMSE <- round(sqrt(mean(unlist(na.omit(et.all))^2)), round.by)
      
      if(RMSE >= RMSE.best){
        no.better <- no.better + 1
      }else{
        no.better <- 0
        RMSE.best <- RMSE
        V.best <- V
        W.best <- W
        C0.best <- C
        mu0.best <- mu
      }
      
      total.steps <- total.steps+1
      Diff.time <- difftime(Sys.time(), start.time, units = 'sec')
      if(Diff.time < 60){
        time.unit <- 'seconds'
      }else{
        Diff.time <- Diff.time/60
        time.unit <- 'minutes'
      }
      print(paste('Total steps:', total.steps, '| Current RMSE:', 
                  round(RMSE,4), '| Best RMSE:', round(RMSE.best,4), 
                  '| Step time:', round(Diff.time,1), time.unit))
      
    }
    
    V <- V.best
    W <- W.best
    C0 <- C0.best
    m0 <- mu0.best
    
    # Save the estimated variance components
    V.name <- paste('V_', stratification.group, sep='')
    V.list.out[[length(V.list.out)+1]] <- V
    names(V.list.out)[length(V.list.out)] <- V.name
    
    W.name <- paste('W_', stratification.group, sep='')
    W.list.out[[length(W.list.out)+1]] <- W
    names(W.list.out)[length(W.list.out)] <- W.name
    
    mu0.name <- paste('mu0_', stratification.group, sep='')
    mu0.list.out[[length(mu0.list.out)+1]] <- mu0
    names(mu0.list.out)[length(mu0.list.out)] <- mu0.name
    
    C0.name <- paste('C0_', stratification.group, sep='')
    C0.list.out[[length(C0.list.out)+1]] <- C0
    names(C0.list.out)[length(C0.list.out)] <- C0.name
    
  }

  return(list(
    'V.list'=V.list.out,
    'W.list' = W.list.out,
    'mu0.list' = mu0.list.out,
    'mu0.list' = C0.list.out,
    'C0.list' = C0.list.out,
    'total.steps' = total.steps
  ))
  
}



## A function for optimizing the discount factor, delta
optimize.delta <- function(deltas=seq(from=0.5, to=1, by=0.01),
                           Data,
                           identifyer,
                           stratify.by,
                           mu0.list, 
                           C0.list,
                           V.list, 
                           relevant.names,
                           Spline.list=NA,
                           time.var){

  # These are the steps we can use to search the parameter space
  steps <- c(0.2, 0.05, 0.01)
  
  # Make an empty data frame, which will be filled up as we iterate over the delta values
  out.all <- data.frame()
  
  # # Calculate the total number of times we need to run the DLM
  # N <- length(deltas)*length(unique(Learningset[,identifyer]))
  # counter <- 0
  
  # If identifyer is NA, we can handle it by adding a dummy column
  if(is.na(identifyer)){
    D.full$Dummy_ID <- 1
    identifyer <- 'Dummy_ID'
  }
  
  # Iterate over the stratification groups
  if(is.na(stratify.by)){
    Data$Dummy <- 1
    stratify.by <- 'Dummy'
  }
  
  Learningset.full <- Data
  best.delta.list <- list()
  for(stratification.group in unique(Learningset.full[,stratify.by])){
    
    Learningset <- subset(Learningset.full, Learningset.full[,stratify.by] == stratification.group)
    
    # Fist we need to see if we have stratification-specific versions of C0
    Grep <- grep(pattern = stratification.group, x = names(C0.list))
    if( length(Grep) == 1 ){
      C0 <- C0.list[[Grep]]
    }else{
      # if not, we just use the first element of C0
      C0 <- C0.list[[1]]
    }

    # Now we need to see if we have stratification-specific versions of mu0
    Grep <- grep(pattern = stratification.group, x = names(mu0.list))
    if( length(Grep) == 1 ){
      mu0 <- mu0.list[[Grep]]
    }else{
      # if not, we just use the first element of C0
      mu0 <- mu0.list[[1]]
    }

    # Fist we need to see if we have stratification-specific versions of C0
    Grep <- grep(pattern = stratification.group, x = names(V.list))
    if( length(Grep) == 1 ){
      V <- V.list[[Grep]]
    }else{
      # if not, we just use the first element of C0
      V <- V.list[[1]]
    }

    for(Step in steps){
      
      # Make smarter delta vectors
      delta_range <- range(deltas)
      deltas <- round( sort(unique(c(delta_range, seq(from=delta_range[1], to=delta_range[2], by=Step)))),2)
      deltas <- c(deltas, deltas[1]+0.01)
      if(nrow(out.all) > 0){
        deltas <- deltas[-which(deltas %in% out.all$delta)]
      }
      deltas <- deltas[-length(deltas)]
      
      
      # Iterate over the delta values
      counter <- 0
      N <- length(deltas)*length(unique(Learningset[,identifyer]))
      for(delta in deltas){
        
        # Define an empty vector to collect all forecast errors in
        et.all <- c()
        
        # Iterate over all unique IDs in the Learningset
        for(ID in unique(Learningset[,identifyer])){
          
          # Extract the relevant subset for the current ID
          ID.set <- subset(Learningset, Learningset[,identifyer] == ID)
          
          # Apply the DLM
          res <- runDLM(Data=ID.set, mu0=mu0, C0=C0, V=V, W=NA, adjust.W=FALSE, 
                        delta=delta, relevant.names=relevant.names, 
                        Spline.list=Spline.list, time.var=time.var)
          
          # Update the counter and report progress (every 5%)
          counter <- counter + 1
          progress <- round(counter/N*100,2)
          
          if (counter %% 100 == 0) {
            print(paste('Optimizing delta, steps:', 
                        Step, 
                        '| Progress:', 
                        progress, '%'))
          }
          
          # Save the forecast errors
          et <- na.omit(unlist(res$et))
          et.all <- c(et.all, et)
        }
        
        # Calculate the root mean squared error for this delta value and save it
        RMSE <- sqrt(mean(et.all^2))
        out.all <- rbind(out.all, cbind(delta, RMSE))
        
      }
      
      # Find the minimum from this run and redefine the delta range for the next Step
      out.all <- out.all[order(out.all$delta),]
      i <- which(out.all$RMSE == min(out.all$RMSE))
      deltas <- out.all$delta[(i-1):(i+1)]
      
    }
    
    
    # Identify the best delta value for this data set
    plot(out.all, type='b')
    best.delta <- out.all$delta[which(out.all$RMSE == min(out.all$RMSE))]
    
    best.delta.list[[length(best.delta.list)+1]] <- best.delta
    
  }

  
  return(best.delta.list)
}


## Extract relevant information from a (multivariate) DLM as a data frame

## Things adapted to this specific model:
## - If the DLM is univariate then the dataframe has to be transposed
## - Unique column names for the multiple cos and sin columns

extract.res <- function(res){
  
  # Get the raw observations
  df_Yt <- data.frame(t(sapply(res$Yt,c)))
  if(nrow(df_Yt) == 1){
      df_Yt <- t(df_Yt)
    }
  colnames(df_Yt) <- relevant.names
  df <- df_Yt
  
  # Get the filtered means 
  df_mt <- data.frame(t(sapply(res$mt,c)))
  if(nrow(df_mt) == 1){
      df_mt <- t(df_mt)
    }
  colnames(df_mt) <- make.unique(paste('Mt', rownames(res$mt[[1]]), sep=''))
  df <- cbind(df, df_mt)
  
  # Get the forecasts
  df_ft <- data.frame(t(sapply(res$ft,c)))
  if(nrow(df_ft) == 1){
      df_ft <- t(df_ft)
    }
  colnames(df_ft) <- paste('Ft', relevant.names,sep='')
  df <- cbind(df, df_ft)
  
  # Get the raw forecasts errors
  df_et <- data.frame(t(sapply(res$et,c)))
  if(nrow(df_et) == 1){
      df_et <- t(df_et)
    }
  colnames(df_et) <- make.unique(paste('Et', rownames(res$et[[1]]),sep=''))
  df <- cbind(df, df_et)
  
  # Get the standardized forecasts errors
  df_ut <- data.frame(t(sapply(res$ut,c)))
  if(nrow(df_ut) == 1){
      df_ut <- t(df_ut)
    }
  colnames(df_ut) <- make.unique(paste('Ut', rownames(res$ut[[1]]),sep=''))
  df <- cbind(df, df_ut)
  
  return(df)
  
}



# ## Validate the DLM by assessing whether the standardized forecast errors follow a standard normal distribution
# validate_DLM <- function(Data, mu0.list, C0.list, V.list, W.list=NA, delta.list=NA, relevant.names, Spline.list, time.var, stratify.by){
# 
#   D.full <- Data
# 
#   # If stratify.by is na, we need to add a Dummy variabe
#   if(is.na(stratify.by)){
#     D$Dummy <- 1
#     stratify.by <- 'Dummy'
#   }
# 
#   # Iterate over the different stratification groups
#   for(stratification.group in unique(D.full[,stratify.by])){
# 
#     # Get the relevant subset of the data
#     D <- subset(D.full, D.full[,stratify.by] == stratification.group)
# 
#     # Fist we need to see if we have stratification-specific versions of C0
#     Grep <- grep(pattern = stratification.group, x = names(C0.list))
#     if( length(Grep) == 1 ){
#       C0 <- C0.list[[Grep]]
#     }else{
#       # if not, we just use the first element of C0.list
#       C0 <- C0.list[[1]]
#     }
# 
#     # Now we need to see if we have stratification-specific versions of mu0
#     Grep <- grep(pattern = stratification.group, x = names(mu0.list))
#     if( length(Grep) == 1 ){
#       mu0 <- mu0.list[[Grep]]
#     }else{
#       # if not, we just use the first element of mu0.list
#       mu0 <- mu0.list[[1]]
#     }
# 
#     # V0 have the option of being NA - we need to handle that, so we arbitrarily use 1/10 of the relevant diagonal values from C0
#     Grep <- grep(pattern = stratification.group, x = names(V.list))
#     if( length(Grep) == 1 ){
#       V <- V.list[[Grep]]
#     }else{
#       # if not, we just use the first element of V.list
#       V <- V.list[[1]]
#     }
# 
#     # W0 have the option of being NA - we need to handle that, so we arbitrarily use 1/10 of the values from C0
#     if(!identical(W.list, NA)){
#       Grep <- grep(pattern = stratification.group, x = names(W.list))
#       if( length(Grep) == 1 ){
#         W <- W.list[[Grep]]
#       }else{
#         # if not, we just use the first element of W.list
#         W <- W.list[[1]]
#       }
#     }
# 
#     # delta have the option of being NA - we need to handle that, so we arbitrarily use 1/10 of the values from C0
#     if(!identical(delta.list, NA)){
#       Grep <- grep(pattern = stratification.group, x = names(delta.list))
#       if( length(Grep) == 1 ){
#         delta <- delta.list[[Grep]]
#       }else{
#         # if not, we just use the first element of W.list
#         delta <- delta.list[[1]]
#       }
#     }
# 
# 
#     #### Check the distribution of the standardized forecast errors on the test set
#     #    - it should be standard normally distributed!
#     ut.all <- cbind()
#     for(ID in unique(D[,identifyer])){
#       ID.set <- subset(Data, D[,identifyer] == ID)
#       
#       print(Spline.list)
# 
#       # We first check to see if we have a W and use if we do - if not, we'll assume we have a discount factor and use that instead
#       if(!identical(W, NA)){
#         res <- runDLM(Data = ID.set, mu0 = mu0,  C0 = C0,  V=V,    W=W,     adjust.W=FALSE, delta=NA,    relevant.names = relevant.names, Spline.list=Spline.list, time.var=time.var, stratify.by)
#       }else{
#         res <- runDLM(Data = ID.set, mu0 = mu0,  C0 = C0,  V=V,    W=NA,    adjust.W=FALSE, delta=delta, relevant.names = relevant.names, Spline.list=Spline.list, time.var=time.var, stratify.by)
#       }
# 
#       # get the standardized forecast errors as a data frame
#       ut <- t(as.data.frame(res$ut))
#       row.names(ut) <- NULL
#       ut.all <- rbind(ut.all, ut)
#     }
#     ut.all <- as.data.frame(ut.all)
#     row.names(ut.all) <- NULL
#     colnames(ut.all) <- relevant.names
# 
#     # - plot it and get the percentage outside of the 95 % CI
#     par(mfrow=c(1,1))
#     for(name in relevant.names){
#       A <- na.omit(ut.all[,name])
#       percent.outside.CI <- round(length(which(A > -1.96 & A < 1.96))/length(A),3) * 100
#       hist(ut.all[,name],breaks = length(A)/10, main=paste(name, '|', percent.outside.CI, '%'), xlim=c(-4,4))
#       abline(v=0, col='red')
#       abline(v=-1.96, col='blue')
#       abline(v=1.96, col='blue')
#     }
# 
#   }
#   return()
# }
#   
# 
# 
# 
# 
