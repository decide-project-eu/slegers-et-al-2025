#### Custom get.Gt function for multivariate DLM with harmonic waves and 
#### trend components:

# NOTICE: all input elements should be vectors, where each element in the 
# vector corresponds to the appropriate value for the corresponding name in 
# relevant.names
# Description of inputs:
# relevant.names is a vector of column names which should be modelled in the 
# DLM
# N.w.all is a vector, where each element is the number of waves used to 
# decribe the corresponding name in relevant.names - the number of waves can 
# also be 0
# trend.all is a vector of boolean values (TRUE or FALSE); if a given value is 
# TRUE, then the corresponding name in relevant.names should be modelled with a 
# trend component
# w.all is a vector of radial frequencies, e.g. (2*pi)/24 for daily patterns 
# measured per hour or (2*pi)/12 for yearly patterns measured per month, etc.; 
# each value describes the period of interest for the corresponding name in 
# relevant.names. Te value can also be NA.

# Examples:
# relevant.names <- c('water', "temp_dry", "humidity", "wind_speed" , "wind_dir")
# N.w.all <- c(3, 1, 1, 1, 0)
# trend.all <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
# w.all <- c( (2*pi)/24,  (2*pi)/24, (2*pi)/24, (2*pi)/24, NA)

get.Gt <- function(N.w.A=N.w.all,
                   trend.A=trend.all,
                   relevant.names.A=relevant.names,
                   w.A=w.all){
  
  # Get all the input variables - as it is now, these are expected to be globally defined. It's not pretty, but it works for now
  w.all=w.A
  relevant.names=relevant.names.A
  N.w.all = N.w.A
  trend.all = trend.A
  
  # Here we determine the number of total rows and columns (and hence the length of the diagonal) in the Gt matrix; 
  # one level per name in relevant names
  # one level per trend
  # two levels per harmonic wave
  N <- length(relevant.names) + length(which(trend.all)) + sum(N.w.all*2)
  
  # Here we make an empty Gt matrix, which is just 1 down the diagonal
  Gt <- diag(1,N)
  
  # j is an index to help us keep track of where the sub-Ft of of each name in relevant name should be placed
  j <- 1
  
  # Now we build the Gt, one name in relevant.names at a time
  for(i in 1:length(relevant.names)){
    
    # Get the values specific to this name
    N.w <- N.w.all[i]
    trend <- trend.all[i]
    w <- w.all[i]
    
    # Add the element for the trend component, if relevant
    if(trend==TRUE){
      j_add <- sum(2+2*N.w)
      Gt[j,j+1] <- 1
      wave.i <- 2
    }else{
      j_add <- sum(1+2*N.w)
      wave.i <- 1
    }
    
    # Add the elements for each wave for this name, if relevant
    if(N.w > 0){
      count <- 0
      wave.elements <- c(cos(w), sin(w), -sin(w), cos(w))
      for(n in 1:N.w){
        Gt[j+wave.i,j+wave.i] <- cos(n*w)
        Gt[j+wave.i,j+wave.i+1] <- sin(n*w)
        Gt[j+wave.i+1,j+wave.i] <- -sin(n*w)
        Gt[j+wave.i+1,j+wave.i+1] <- cos(n*w)
        wave.i <- wave.i + 2
      }
    }
    
    # Update j
    j <- j + j_add
    
  }
  
  
  return(Gt)
}


# Custom get.Ft function for the "harmonic waves" method
get.Ft <- function(N.w.A=N.w.all,
                   trend.A=trend.all,
                   relevant.names.A=relevant.names,
                   w.A=w.all){
  
  # Get all the input variables - as it is now, these are expected to be globally defined. It's not pretty, but it works for now
  w.all=w.A
  relevant.names=relevant.names.A
  N.w.all = N.w.A
  trend.all = trend.A
  
  # Here we determine the number of total rows in the Ft matrix; 
  # one level per name in relevant names
  # one level per trend
  # two levels per harmonic wave
  N <- length(relevant.names) + length(which(trend.all)) + sum(N.w.all*2)
  
  # Ft.base is a vector of zeros, to which the proper values will be added for each name in relevant.names
  Ft.base <- rep(0,N)
  
  # j is an index to help us keep track of where the sub-Ft of of each name in relevant name should be placed
  j <- 1
  
  # Now we build the Ft, one name in relevant.names at a time
  Ft <- cbind()
  for(i in 1:length(relevant.names)){

    # Get the values specific to this name
    N.w <- N.w.all[i]
    trend <- trend.all[i]
    w <- w.all[i]
    
    # Add the element for the trend component, if relevant
    if(trend == TRUE){
      Ft.i <- c(1,0)
    }else{
      Ft.i <- c(1)
    }
    
    # Add the elements for each wave for this name, if relevant
    if(N.w > 0){
      for(n in 1:N.w){
        Ft.i <- c(Ft.i, c(1,0))
      }
    }
    
    # Add the sub-Ft for this name to the blank Ft.base and update the final Ft
    Ft.j <- Ft.base
    Ft.j[j:(j+length(Ft.i)-1)] <- Ft.i
    Ft <- rbind(Ft, Ft.j)
    
    # Update j
    j <- j + length(Ft.i)
    
  }
  
  # Turn Ft into a matrix
  Ft <- t(as.matrix(Ft))
  
  return(Ft)
  
}



# Function for defining the LM for a data set as the sum of some trend and a number of harmonic waves
make.LM.wHarmonics <- function(Data, w, N.w, trend, relevant.name, time.var.name, plot.it=FALSE, remove.zeros=FALSE, round.by=3){

  # Remove zero-values, if relevant
  if(remove.zeros == TRUE){
    remove.i <- which(D[,relevant.name] == 0)
    D <- D[-remove.i,]
  }
  
  # We make an empty vector, into which we will add the relevant text for the linear function
  lm.vector <- c()
  
  # We also make an empty vector, into which we will add the names for the parameter vector
  names.vector <- c()
  
  # We add the target variable (the variable we model with this collection of harmonics)
  lm.vector <- c(lm.vector, relevant.name, '~')
  names.vector <- c(names.vector, relevant.name)
  
  # Make the linear trend component, i relevant
  if(trend == TRUE){
    lm.vector <- c(lm.vector, time.var.name)
    names.vector <- c(names.vector, 'Trend')
  }
  
  # Now we add the relevant number of harmonics, one at a time
  if(N.w > 0){
    for(i in 1:N.w){
      a <- paste('+ cos(', i, '*w*', time.var.name, ') + sin(', i, '*w*', time.var.name, ')')
      lm.vector <- c(lm.vector, a)
      names.vector <- c(names.vector, paste('wave.', i, '.cos', sep=''), paste('wave.', i, '.sin', sep=''))
    }
  }
  
  # We can now make the final name vector
  names.vector <- paste(names.vector, sep='')
  
  # If there is no trend and no waves, we need to add a NULL to the model
  if(length(lm.vector)==2){
    lm.vector <- c(lm.vector, 'NULL')
  }

  # We only consider the data in the beginning
  # D.A <- subset(D, D[,time.var.name] <= 0.1*max(D[,time.var.name]))
  
  # Now we make it into a linear function - without random effect
  b1 <- as.formula(paste(lm.vector, collapse=' '))
  lm.1 <- lm(b1, data = Data )
  S1 <- summary(lm.1)

  # - now we directly have the initial parameter vector
  # mu0 <- matrix(S1$coefficients[,'Estimate'])
  mu0 <- matrix(S1$coefficients[,'Estimate'])
  mu0 <- matrix(mu0[1:length(names.vector)])
  rownames(mu0) <- names.vector
  
  # - and we directly have the initial prior variance matrix
  # C0 <- as.matrix(vcov(lm.1))
  C0 <- as.matrix(vcov(lm.1))
  C0 <- C0[1:length(names.vector),1:length(names.vector)]
  
  # Get the adjusted R^2 for the linear model
  S <- summary(lm.1)
  adj.r.squared <- S$adj.r.squared

  
  # Return the mu0 and C0 for use in the DLM
  return(list('mu0'=mu0,
              'C0'=C0,
              'LM'=lm.1,
              'adj.r.squared'=adj.r.squared))
  
}



# relevant.names <- c("temp_dry", "humidity", "wind_speed" , "wind_dir")
# N.w.all <- c(1, 1, 1, 0)
# trend.all <- c(FALSE, FALSE, FALSE, FALSE)
# w.all <- c( (2*pi)/24, (2*pi)/24, (2*pi)/24, NA)
# time.var.names.All <- c("Hour", "Hour", "Hour", "Hour")
get.mu0_and_C0 <- function(Data, relevant.names, time.var.names.All, trend.all, N.w.all, w.all){
  
  N <- length(relevant.names) + length(which(trend.all)) + sum(N.w.all*2)
  mu0.all <- rep(0, N)
  C0.all <- diag(0, N)
  
  j <- 1
  names.all <- c()
  for(i in 1:length(relevant.names)){
    
    # Get the values specific to this name
    relevant.name <- relevant.names[i]
    N.w <- N.w.all[i]
    trend <- trend.all[i]
    w <- w.all[i]
    time.var.name <- time.var.names.All[i]
    
    a <- make.LM.wHarmonics(Data, w, N.w, trend, relevant.name, time.var.name, plot.it=FALSE, remove.zeros=FALSE, round.by=3)
    mu0 <- a$mu0
    C0 <- a$C0
    
    mu0.all[j:(j+length(mu0)-1)] <- mu0
    C0.all[j:(j+length(mu0)-1),j:(j+length(mu0)-1)] <- C0
    
    names.all <- c(names.all, rownames(mu0))
    
    j <- j + length(mu0)
    
  }
  
  mu0.all <- matrix(mu0.all)
  rownames(mu0.all) <- names.all
  
  C0.all <- as.matrix(C0.all)
  rownames(C0.all) <- names.all
  colnames(C0.all) <- names.all
  
  return(list('mu0'=mu0.all,
              'C0'=C0.all))
}












