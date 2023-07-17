convert<- function(WD){
  ################################################################
  ## Maps the angles from [-pi, pi] to the bearing [0,2*pi]
  ################################################################
  # Map the angles to [0,2pi]
  WD<- ifelse(WD < 0 ,2 * pi + WD, WD)
  # Change the angles to the navigational-bearing scale
  return(WD.bearing<- (pi/2 - WD) %% (2 * pi) )
  #return(WD.bearing<- (WD)%%(2*pi) )
  #data$WD.bearing<- WD.bearing
  #return(WD)
}

bin<- function(WD, type, width, n.knots){
  ############################################################
  ## Bins the data with respect to wind direction
  ############################################################
  ## WD = wind direction data
  ## type = ew.bin/ef.bin choice btw equal width and equal freq. bins
  ## width = width of the bin in terms of degrees for ew.bin
  ## n.knots = number of knots for ef.bins
  if(type== "ew.bin"){
    # Length of a bin
    deg<- width*pi/180
    # Define the range of each bin
    dir.breaks <- seq(0,  2*pi, deg)
    # Assign each direction to a bin range
    dir.binned <- cut(WD, breaks = dir.breaks, ordered_result = TRUE)
    # Generate labels
    dir.labels <- as.character(c(seq(0, 2*pi-deg, by = deg), 0))
    # replace ranges with bin labels
    levels(dir.binned) <- dir.labels
    # Assign bin names to the original data set
  return(BIN <- dir.binned)
  }
  if(type == "ef.bin"){
    # Define the knots
    degf <- n.knots
    knots<- quantile(WD, 1:degf/(degf+1))
    # Define the range of each bin
    dir.breaks <- c(0, knots, 2*pi)
    # Assign each direction to a bin range
    dir.binned <- cut(WD, breaks = dir.breaks, ordered_result = TRUE)
    # Generate labels
    dir.labels <- as.character(c(0, knots, 0))
    # replace ranges with bin labels
    levels(dir.binned) <- dir.labels
    # Assign bin names to the original data set
    return(BIN <- dir.binned)
  }
}

Harmonic <- function(theta, K){
############################################################
#### Fourier series
############################################################
  t <- outer(theta, 1:K)
  return(cbind(apply(t, 2, cos), apply(t, 2, sin)))
}

X.median.func = function(WD, K, type, width, n.knots){
############################################################
#### Computes the sine/cosine of the wind direction median
############################################################
  # WD = wind direction data
  # K = number of pairs of Fourier series
  # type = ew.bin for equal width bin
  # type = ef.bins for equal frequency bins
  # width = the width (in degrees) of each bin, when equal width bins are chosen
  # n.knots = number of knots when equal frequency is chosen
  bin.wd = bin(WD = WD, type = type, width = width, n.knots = n.knots) 
  df = data.frame(WD = WD, Bins = bin.wd)
  #Compute the Fourier series
  dir.median = aggregate(df$WD, by = list(df$Bins), function(x) summary(x))$x[,3]
  return(X.median = Harmonic(dir.median, K ))
}

binned.pp.step1 = function(WD, WS, type, width, n.knots, th ){
############################################################
#### First step of the algorithm using the point process: binning of WD
############################################################
  # WD = wind direction data
  # WS = wind speed data
  # type = ew.bin for equal width bin
  # type = ef.bins for equal frequency bins
  # width = the width (in degrees) of each bin, when equal width bins are chosen
  # n.knots = number of knots when equal frequency is chosen
  # th = threshold to be aplied in each bin

  # Bin WD
  bin.wd = bin(WD = WD, type = type, width = width, n.knots = n.knots) 
  # Organize data in a dataframe
  df = data.frame(WS = WS, WD = WD, Bins = bin.wd)
  ## In each bin fit a point process
  ppt.fit = aggregate(df$WS, by = list(df$Bins), 
                      function(x) pp.fit(x, threshold = quantile(x, th),npy = 92, show = F))
  # Save the parameter estimates and their se
  loc.est = as.numeric(sapply(ppt.fit$x[,12], function(x) x[1])) 
  scale.est = as.numeric(sapply(ppt.fit$x[,12], function(x) x[2])) 
  shape.est = as.numeric(sapply(ppt.fit$x[,12], function(x) x[3])) 
  thresh.est = as.numeric(sapply(ppt.fit$x[,4], function(x) x[1])) 
  se.loc.est = as.numeric(sapply(ppt.fit$x[,14], function(x) x[1])) 
  se.scale.est = as.numeric(sapply(ppt.fit$x[,14], function(x) x[2])) 
  se.shape.est = as.numeric(sapply(ppt.fit$x[,14], function(x) x[3])) 
  binned.param = list(cbind(loc.est, scale.est, shape.est, thresh.est), cbind(se.loc.est, se.scale.est, se.shape.est))
  return(binned.param)
}
                                   
binned.pp.step2 = function(binned.param, X.median,K, RL = c(5, 20, 50)){
############################################################
#### Second step of the algorithm: dependence of the point ptocess' paramters on WD
############################################################
  # binned.param = parameters from step 1 as a list of two lists, where 1st list
  #contains the mle's and second list contains the se's
  # X.median = Fourier series of the binned wind directions medians
  # K = number of pairs of Fourier series to be used
  # RL = return levels to be calculated
  
  ## Regress the parameter estimates on WD
  loc.param = lm(binned.param[[1]][,1] ~ X.median, weights = 1/(binned.param[[2]][,1])^2)
  scale.param = lm(log(binned.param[[1]][,2]) ~ X.median, weights = 1/(log(binned.param[[2]][,2]))^2)
  shape.param = lm(binned.param[[1]][,3] ~ X.median, weights = 1/(binned.param[[2]][,3])^2)
  thresh.param = lm(binned.param[[1]][,4] ~ X.median)#, weights = 1/(binned.param[[1]][[2]][,3])^2)
  
  # Compute the parameters of the point process
  xg = seq(0, 2*pi, 0.01)
  yg.loc.param = cbind(rep(1, length(xg)), Harmonic(xg, K)) %*% loc.param$coefficients
  yg.scale.param = cbind(rep(1, length(xg)), Harmonic(xg, K)) %*% scale.param$coefficients
  yg.shape.param = cbind(rep(1, length(xg)), Harmonic(xg, K)) %*% shape.param$coefficients
  yg.thresh.param = cbind(rep(1, length(xg)), Harmonic(xg, K)) %*% thresh.param$coefficients
  
  ## Compute the RLs
  RL.pp.est = array(dim = c(629, 3))
  for(i in 1:629){ 
    #RL.pp.est[i, ] = ppq(c(yg.loc.param[i], yg.scale.param[i], yg.shape.param[i]), 
    #yg.thresh.param[i], 92,  1 / (RL * 92 * (1 - 0.95)))
    RL.pp.est[i, ] = qevd( 1 - 1 / (RL * 92 * (1 - 0.95)), loc = yg.loc.param[i],
                           scale = exp(yg.scale.param[i]), shape = yg.shape.param[i],
                           threshold = yg.thresh.param[i], type = "PP")
  }
  return(RL.pp.est)
}

binned.gpd.step1 = function(WD, WS, type, width, n.knots, th ){
############################################################
#### First step of the algorithm using the GPD: binning of WD
############################################################
  ### arguments are the same as for point process
  # Bin WD
  # type = ew.bin for equal width bin
  # type = ef.bins for equal frequency bins
  bin.wd = bin(WD = WD, type = type, width = width, n.knots = n.knots) 
  df = data.frame(WS = WS, WD = WD, Bins = bin.wd)
  ## In each bin fit a point process
  gpd.fit = aggregate(df$WS, by = list(df$Bins), 
                      function(x) gpd.fit(x, threshold = quantile(x, th), 
                                          npy = 92, show = F))
  # Save the parameter estimates and their se
  scale.est = as.numeric(sapply(gpd.fit$x[,10], function(x) x[1])) 
  shape.est = as.numeric(sapply(gpd.fit$x[,10], function(x) x[2])) 
  thresh.est = as.numeric(sapply(gpd.fit$x[,4], function(x) x[1])) 
  se.scale.est = as.numeric(sapply(gpd.fit$x[,13], function(x) x[1])) 
  se.shape.est = as.numeric(sapply(gpd.fit$x[,13], function(x) x[2])) 
  binned.param = list(cbind(scale.est, shape.est, thresh.est), cbind(se.scale.est, se.shape.est))
  return(binned.param)
}

binned.gpd.step2 = function(binned.param, X.median, K, RL = c(5, 20, 50)){
############################################################
#### Second step of the algorithm: dependence of the GPD's paramters on WD
############################################################
  ### arguments are the same as for point process
  
  ## Regress the parameter estimates on WD
  scale.param = lm(log(binned.param[[1]][,1]) ~ X.median, weights = 1/(log(binned.param[[2]][,1]))^2)
  shape.param = lm(binned.param[[1]][,2] ~ X.median, weights = 1/(binned.param[[2]][,2])^2)
  thresh.param = lm(binned.param[[1]][,3] ~ X.median)#, weights = 1/(binned.param[[1]][[2]][,3])^2)
  
  # Compute the parameters of the point process
  xg = seq(0, 2*pi, 0.01)
  yg.scale.param = cbind(rep(1, length(xg)), Harmonic(xg, K)) %*% scale.param$coefficients
  yg.shape.param = cbind(rep(1, length(xg)), Harmonic(xg, K)) %*% shape.param$coefficients
  yg.thresh.param = cbind(rep(1, length(xg)), Harmonic(xg, K)) %*% thresh.param$coefficients
  
  ## Compute the RLs
  RL.gpd.est = array(dim = c(629, 3))
  for(i in 1:629){ 
    RL.gpd.est[i, ] = gpdq(c(exp(yg.scale.param[i]), yg.shape.param[i]), 
    yg.thresh.param[i], 1 / (RL * 92 * (1 - 0.95)))
    #RL.gpd.est[i, ] = qevd(1 - 1 / (RL * 92 * (1 - 0.95)), scale = exp(yg.scale.param[i]), 
                           #shape = yg.shape.param[i],threshold = yg.thresh.param[i],
                          # type = "GP")
  }
  return(RL.gpd.est)
}

binned.gev.step1 = function(WD, WS, type, width, n.knots){
############################################################
#### First step of the algorithm using the GEV distribution: binning of WD
############################################################
  ### arguments are the same as for point process
  
  # In each bin find the maximum
  nyears = 55
  WS.binned.max = bin.wd = list()
  for(i in 1:nyears){
    bin.wd[[i]] = bin(WD = WD[,i], type = type, width = width, n.knots = n.knots) 
    df = data.frame(WS = WS[,i], WD = WD[,i], Bins =  bin.wd[[i]])
    WS.binned.max[[i]] = aggregate(df$WS, by = list(df$Bins), max)$x
  }
  list.length = array()
  for(i in 1:nyears){
    list.length[i] = length(WS.binned.max[[i]])
  }
  ## Fit a GEV for the binned maxs
  gev1.fit = list()
  for(i in 1:max(list.length)){
    d = sapply(WS.binned.max, function(x) x[i])
    d = d[!is.na(d)]
    gev1.fit[[i]] =  gev.fit(d, show = F)
  }
  loc.mle = loc.mle.se = scale.mle = scale.mle.se = shape.mle = shape.mle.se = array()
  for(i in 1:max(list.length)){
    loc.mle[i] = gev1.fit[[i]]$mle[1]
    loc.mle.se[i] =gev1.fit[[i]]$se[1]
    scale.mle[i] = gev1.fit[[i]]$mle[2]
    scale.mle.se[i] = gev1.fit[[i]]$se[2]
    shape.mle[i] = gev1.fit[[i]]$mle[3]
    shape.mle.se[i] = gev1.fit[[i]]$se[3]
  }
  return(list(cbind(loc.mle, scale.mle, shape.mle), cbind(loc.mle.se, scale.mle.se, shape.mle.se)))
}

binned.gev.step2 = function(binned.param, X.median, K, RL = c(5, 20, 50)){
############################################################
#### Second step of the algorithm: dependence of the GEV distribution's paramters on WD
############################################################
### arguments are the same as for point process
  
  loc.param.pre.hc = lm(binned.param[[1]][,1] ~ X.median, weights = 1/(binned.param[[2]][,1])^2)
  scale.param.pre.hc = lm(log(binned.param[[1]][,2]) ~ X.median, weights = 1/(log(binned.param[[2]][,2]))^2)
  shape.param.pre.hc = lm(binned.param[[1]][,3] ~ X.median, weights = 1/(binned.param[[2]][,3])^2)
  xg = seq(0, 2*pi, len = 629)
  yg.loc.hc <- cbind(rep(1, 629), Harmonic(xg, K)) %*% loc.param.pre.hc$coefficients
  yg.scale.hc <- cbind(rep(1, 629), Harmonic(xg, K)) %*% scale.param.pre.hc$coefficients
  yg.shape.hc <- cbind(rep(1, 629), Harmonic(xg, K)) %*% shape.param.pre.hc$coefficients
  
  RL.gev.est = array(dim = c(629, 3))
  for(i in 1:629) {
    #RL.gev.est[i,] = qevd(1 - 1 / RL, loc = yg.loc.hc[i], 
                         # scale = exp(yg.scale.hc[i]), 
                          #shape = yg.shape.hc[i], type = "GEV")
    RL.gev.est[i,] = gevq(c(yg.loc.hc[i], 
                           exp(yg.scale.hc[i]), 
                          yg.shape.hc[i]), 1 / RL)
  }
  return(RL.gev.est)
}
