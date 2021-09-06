



# 1.1 Function of True (Naive) Model  --- Yule-Walker -------------------------------

ARyw <- function(ts, order) {
  # xm <- colMeans(ts)
  # x <- sweep(x, 2L, xm, check.margin=FALSE)
  GAMMA <- NULL
  xacf <- acf(ts, type = "covariance", lag.max = order, plot = FALSE,
              demean = T)$acf
  gammahat <- xacf[-1]
  gammarow <- xacf[-length(xacf)]
  xacfcan <- c(rev(gammarow),gammarow[-1])
  nGamma <- length(xacf)-1
  for (i in nGamma:1){
    GAMMA <- rbind(GAMMA,xacfcan[i:(i+nGamma-1)])
  }
  phihat <- solve(GAMMA,tol=1e-250) %*% gammahat
  sigma_squ <- gammarow[1] - 2 * sum( phihat *  gammahat) + t(phihat) %*% GAMMA %*% as.matrix(phihat)
  phi0 <- (1-sum(phihat))* mean(ts)
  return(list(phi0=phi0,phihat=phihat,sigma_squ=sigma_squ))
}


ARywClassic <- function(ts, order, alpha0, alpha1, sigma_e){
  GAMMA <- NULL
  xacf <- acf(ts, type = "covariance", lag.max = order, plot = FALSE,
              demean = T)$acf
  gammahat <- xacf[-1]
  gammarow <- xacf[-length(xacf)]
  gammarow[1] <- gammarow[1] - sigma_e^2
  gammarow <- gammarow
  xacfcan <- c(rev(gammarow),gammarow[-1])
  nGamma <- length(xacf)-1
  for (i in nGamma:1){
    GAMMA <- rbind(GAMMA,xacfcan[i:(i+nGamma-1)])
  }
  phihat <- solve(GAMMA,tol=1e-250) %*% gammahat
  sigma_squ <- gammarow[1] - 2 * sum( phihat *  gammahat) + t(phihat) %*% GAMMA %*% as.matrix(phihat)
  phi0 <- (1-sum(phihat))* (mean(ts) - alpha0) / (alpha1)
  return(list(phi0=phi0,phihat=phihat,sigma_squ=sigma_squ))
}

ARywMultiplicative <- function(ts, order, beta0, sigma_u){
  mu <- mean(ts)
  GAMMA <- NULL
  xacf <- acf(ts, type = "covariance", lag.max = order, plot = FALSE,
              demean = T)$acf
  gammahat <- xacf[-1]/beta0^2
  gammarow <- xacf[-length(xacf)]
  gammarow <- gammarow/beta0^2
  gammarow[1] <- gammarow[1]/(sigma_u^2+1) - sigma_u^2/(sigma_u^2+1)*(mu^2)
  
  xacfcan <- c(rev(gammarow),gammarow[-1])
  nGamma <- length(xacf)-1
  for (i in nGamma:1){
    GAMMA <- rbind(GAMMA,xacfcan[i:(i+nGamma-1)])
  }
  phihat <- solve(GAMMA,tol=1e-250) %*% gammahat
  sigma_squ <- gammarow[1] - 2 * sum( phihat *  gammahat) + t(phihat) %*% GAMMA %*% as.matrix(phihat)
  phi0 <- (1-sum(phihat))* mu / beta0 
  return(list(phi0=phi0,phihat=phihat,sigma_squ=sigma_squ))
}

###
# timeseries <- log(Hubei_confirm$cases+0.0001)
# object <- ARyw(log(Hubei_confirm$cases+0.0001), 2)
# n.ahead <- 5


tsbootstrap <- function(ts, FUN, seed=1, Nblock, nBoot=500, order, 
                        alpha0 = NULL, alpha1 = NULL, beta0 = NULL, sigma_e = NULL, sigma_u = NULL ) {
  phi0boot <- NULL
  phihatboot <- NULL
  sigma_squboot <- NULL
  
  Tlength <- length(ts)
  # if (Tlength %% Nblock !=0 ) {stop("The time series cannot break into equal length, try different Nblock")}
  
  blength <- ceiling(Tlength/Nblock)
  
  poollength <- Tlength - blength +1
  
  set.seed(seed)
  
  if (FUN=="ARyw") { arg <- list(ts=ts,order=order)}
  if (FUN=="ARywClassic") { arg <- list(ts=ts,order=order, alpha0=alpha0, alpha1=alpha1, sigma_e=sigma_e)}
  if (FUN=="ARywMultiplicative") { arg <- list(ts=ts,order=order, beta0=beta0, sigma_u=sigma_u)}
  
  for (i in 1:nBoot) {
    indexhead <- sample(poollength, size = Nblock, replace = T)
    tsboot <- unlist(lapply(indexhead, FUN=function(x) { return(ts[x:(x+blength-1)])}))
    tsboot <- tsboot[1:length(ts)]
    
    if (FUN=="ARyw") { arg <- list(ts=tsboot,order=order)} else{
      if (FUN=="ARywClassic") { arg <- list(ts=tsboot,order=order, alpha0=alpha0, alpha1=alpha1, sigma_e=sigma_e)} else {
        if (FUN=="ARywMultiplicative") { arg <- list(ts=tsboot,order=order, beta0=beta0, sigma_u=sigma_u)}
      }
    }
    
    results <- do.call(FUN, arg=arg)
    phi0boot <- c(phi0boot,results$phi0)
    phihatboot <- rbind(phihatboot,drop(results$phihat))
    sigma_squboot <- c(sigma_squboot,results$sigma_squ)
  }
  

  
  return(list(phi0boot.sd=sd(phi0boot), phihatboot.sd= apply(as.matrix(phihatboot), MARGIN = 2, FUN = sd), sigma_squ.sd= sd(sigma_squboot)))  
}





prediction.AR <- function(timeserires, n.ahead, object, order = 1){
  tscandidate <- timeserires
  if (order==1) {
    for (i in 1:n.ahead){
      new.pred <- object$phi0 +  object$phihat * tscandidate[length(tscandidate)]
      tscandidate <- c(tscandidate, new.pred)
    }
  } else {
    for (i in 1:n.ahead){
      new.pred <- object$phi0 +  tscandidate[(length(tscandidate)-order+1):length(tscandidate)] %*% as.matrix(object$phihat)
      tscandidate <- c(tscandidate, new.pred)
    }
  }
  
  return(tscandidate)
}


FittedValue.AR <- function(timeserires, diff, object, order = 1) {
  fittedvalue <- NULL
  if (diff == F) {
    for (t in (order+1):length(timeserires)) {
      pred1 <- prediction.AR(timeserires[(t-order):(t-1)], 1, object, order =order)
      fittedvalue <- c(fittedvalue, pred1[(order+1):(length(pred1))])
    }
    return(fittedvalue)
  } else {
    timeseriresdiff <- diff(timeserires)
    for (t in (order+1):length(timeseriresdiff)) {
      pred1 <- prediction.AR(timeseriresdiff[(t-order):(t-1)], 1,  object, order = order)
      fittedvalue <- c(fittedvalue, pred1[(order+1):(length(pred1))])
    }
    return(fittedvalue + timeserires[(order+1):(length(timeserires)-1)])
  }
}
# EstPred.Case2.1_Fitted <- FittedValue.AR(FataRate_Case3_ON_comp[(Case1tday-laghere):(length(FataRate_Case3_ON_comp)-5)] * 0.54, 
#                                          object = model2.1, diff = F, order = laghere)

# EstPred.Case2.1_Fitted <- FittedValue.AR(FataRate_Case2_ON_comp[(Case1tday-laghere-1):(length(FataRate_Case2_ON_comp)-5)] * 0.54, 
#                                          object = model2.1, diff = T, order = laghere)
# FittedValue.AR(FataRate_Case1_BC_comp[(dateCase1t-laghere-1):length(FataRate_Case1_BC_comp)] * 0.54, object = model2.1, diff = T, order = 1)
# FittedValue.AR(FataRate_Case1_ON_comp[(Case1tday-laghere-1):(length(FataRate_Case1_ON_comp)-5)] * 0.54,
#                object = model2.1, diff = T, order = 4)


pridiction.Kalman <- function(timeseries,n.ahead,object){
  model <- makeARIMA(phi = object$phihat, theta = numeric(0), Delta = numeric(0))
  kf <- KalmanRun(y = ts(timeseries), mod = model, update =T) 
  mod <- attr(kf, "mod")
  
  rsd <- kf$resid - object$phi0
  xtsp <- tsp(ts(rsd))
  n <- length(rsd)
  
  xm <- rep(1, n.ahead) * object$phi0
  
  z <- KalmanForecast(n.ahead, mod)
  pred <- ts(z[[1L]] + xm, start = xtsp[2L] + deltat(rsd), 
             frequency = xtsp[3L])
  se <- ts(sqrt(z[[2L]] * drop(object$sigma_e_squ)), start = xtsp[2L] + 
               deltat(rsd), frequency = xtsp[3L])
  list(pred = pred, se = se)
}


getActive <- function(confirm, death, recovery, I0) {
  ActiveCases <- cumsum(confirm) -  cumsum(death) -  (recovery) + I0
  return(ActiveCases)
}

predictionIntAR1 <- function(h, model) {
  psi <- 1 
  for (i in 1:h){
    psi <- c(psi, model$phihat^i)
  }
  
  Sigmah <- drop(model$sigma_squ) * cumsum(psi^2) 
  bound <- 1.96 * sqrt(Sigmah[1:h])
  return(c(0,bound))
}

predictionIntAR1_add <- function(h, model,sigma_e,alpha1) {
  psi <- 1 
  for (i in 1:(h-1)){
    psi <- c(psi, model$phihat^i)
  }
  
  add <- (drop(model$phiha))^(c(1:h*2)) * sigma_e^2 / alpha1^2
  
  Sigmah <- drop(model$sigma_squ) * cumsum(psi^2) + add
  bound <- 1.96 * sqrt(Sigmah)
  return(c(0,bound))
}


predictionIntAR1_mult <- function(h, model, sigma_u, mu) {
  psi <- 1 
  for (i in 1:(h-1)){
    psi <- c(psi, model$phihat^i)
  }
  
  add <- drop(model$phihat)^(1:h*2) * sigma_u^2 * ( drop(model$sigma_squ)/(1-drop(model$phihat)) + mu^2)
  
  Sigmah <- drop(model$sigma_squ) * cumsum(psi^2) + add
  bound <- 1.96 * sqrt(Sigmah)
  return(c(0,bound))
}

predictionIntAR2 <- function(h, model) {
  psi <- c(1, model$phihat[1])
  if (h>2){
    for (i in 3:h){
      psi <- c(psi,model$phihat[1]*psi[length(psi)] + model$phihat[2] * psi[length(psi) - 1] )
    }
  }
  
  
  Sigmah <- drop(model$sigma_squ) * cumsum(psi^2) 
  bound <- 1.96 * sqrt(Sigmah[1:h])
  return(c(0,bound))
}

predictionIntAR3 <- function(h, model) {
  psi <- c(1, model$phihat[1], model$phihat[1]^2 * model$phihat[2])
  if (h>3){
    for (i in 4:h){
      psi <- c(psi,model$phihat[1]*psi[length(psi)] + model$phihat[2] * psi[length(psi) - 1] + model$phihat[3] * psi[length(psi) - 2] )
    }
  }
  
  
  Sigmah <- drop(model$sigma_squ) * cumsum(psi^2) 
  bound <- 1.96 * sqrt(Sigmah[1:h])
  return(c(0,bound))
}

predictionIntAR2_add <- function(h, model,sigma_e,alpha1) {
  psi <- c(1, model$phihat[1])
  if (h>2){
    for (i in 3:h){
      psi <- c(psi,model$phihat[1]*psi[length(psi)] + model$phihat[2] * psi[length(psi) - 1] )
    }
  }

  add <- c(sum((model$phihat)^2), 
           model$phihat[1]^2*sum(model$phihat^2) + model$phihat[2]^2) * sigma_e^2 / alpha1^2
  
  for (i in 3:h){
    add <- c(add, model$phihat[1]^2 * add[length(add)] + model$phihat[2]^2 * add[length(add)-1] )
  }
  
  Sigmah <- drop(model$sigma_squ) * cumsum(psi^2) + add
  bound <- 1.96 * sqrt(Sigmah)
  return(c(0,bound))
}


predictionIntAR2_mult <- function(h, model, sigma_u, mu) {
  psi <- c(1, model$phihat[1])
  if (h>2){
    for (i in 3:h){
      psi <- c(psi,model$phihat[1]*psi[length(psi)] + model$phihat[2] * psi[length(psi) - 1] )
    }
  }

  add <- c(sum((model$phihat)^2), 
           model$phihat[1]^2*sum(model$phihat^2) + model$phihat[2]^2) * 
    sigma_u^2 * ( drop(model$sigma_squ)/(1-sum(drop(model$phihat))) + mu^2)
  
  for (i in 3:h){
    add <- c(add, model$phihat[1]^2 * add[length(add)] + model$phihat[2]^2 * add[length(add)-1] )
  }
  
  Sigmah <- drop(model$sigma_squ) * cumsum(psi^2) + add
  bound <- 1.96 * sqrt(Sigmah)
  return(c(0,bound))
}

predictionIntAR3_add <- function(h, model,sigma_e,alpha1) {
  psi <- c(1, model$phihat[1], model$phihat[1]^2 * model$phihat[2])
  if (h>3){
    for (i in 4:h){
      psi <- c(psi,model$phihat[1]*psi[length(psi)] + model$phihat[2] * psi[length(psi) - 1] + model$phihat[3] * psi[length(psi) - 2] )
    }
  }
  
  add <- c(sum((model$phihat)^2), 
           model$phihat[1]^2*sum(model$phihat^2) + (model$phihat[2]^2+model$phihat[3]^2),
           model$phihat[1]^2*(model$phihat[1]^2*sum(model$phihat^2) + (model$phihat[2]^2+model$phihat[3]^2)) + model$phihat[2]^2*sum(model$phihat^2) + model$phihat[3]^2
           ) * sigma_e^2 / alpha1^2
  
  for (i in 4:h){
    add <- c(add, model$phihat[1]^2 * add[length(add)] + model$phihat[2]^2 * add[length(add)-1] + model$phihat[3]^2 * add[length(add)-2])
  }
  
  Sigmah <- drop(model$sigma_squ) * cumsum(psi^2) + add
  bound <- 1.96 * sqrt(Sigmah)
  return(c(0,bound))
}


predictionIntAR3_mult <- function(h, model, sigma_u, mu) {
  psi <- c(1, model$phihat[1], model$phihat[1]^2 * model$phihat[2])
  if (h>3){
    for (i in 4:h){
      psi <- c(psi,model$phihat[1]*psi[length(psi)] + model$phihat[2] * psi[length(psi) - 1] + model$phihat[3] * psi[length(psi) - 2] )
    }
  }
  
  add <- c(sum((model$phihat)^2), 
           model$phihat[1]^2*sum(model$phihat^2) + (model$phihat[2]^2+model$phihat[3]^2),
           model$phihat[1]^2*(model$phihat[1]^2*sum(model$phihat^2) + (model$phihat[2]^2+model$phihat[3]^2)) + model$phihat[2]^2*sum(model$phihat^2) + model$phihat[3]^2
          ) * sigma_u^2 * ( drop(model$sigma_squ)/(1-sum(drop(model$phihat))) + mu^2)
  
  for (i in 4:h){
    add <- c(add, model$phihat[1]^2 * add[length(add)] + model$phihat[2]^2 * add[length(add)-1] + model$phihat[3]^2 * add[length(add)-2])
  }
  
  Sigmah <- drop(model$sigma_squ) * cumsum(psi^2) + add
  bound <- 1.96 * sqrt(Sigmah)
  return(c(0,bound))
}



## Auxiliary functions
imputezero <- function(tseries){
  tseries <- unlist(tseries)
  zerostrike <- rep(NA, length(tseries)+1)
  zerostrike[length(tseries)+1] <- 0
  for (i in length(tseries):1){
    if (tseries[i] == 0) {zerostrike[i] = zerostrike[i+1] + 1} else zerostrike[i] = 0
  }
  zerostrike[zerostrike>3] = 0
  for (i in 2:length(tseries)){
    if (zerostrike[i] > zerostrike[i-1]) {
      tseries[(i-1):(i-1+zerostrike[i])] = round(tseries[i-1]/(zerostrike[i]+1))}
  }
  return(tseries)
} 




nodiff_fit <- function (tseries, lag = 1, sigma_e, sigma_u) {
  model1 <- ARyw(tseries, order = lag)
  model2.1 <- ARywClassic(tseries, order = lag, alpha1 =  1/0.54, alpha0 = 0, sigma_e = sigma_e) 
  model2.2 <- ARywClassic(tseries, order = lag, alpha1 =  1/0.54, alpha0 = 0, sigma_e = sigma_e*2) 
  
  model3.1 <- ARywMultiplicative(tseries, order = lag, beta0 = 1/0.54, sigma_u = sigma_u)  
  model3.2 <- ARywMultiplicative(tseries, order = lag, beta0 = 1/0.54, sigma_u = sigma_u*2) 
  
  
  sdmodel1   <- tsbootstrap(FUN = "ARyw", ts = tseries, order = lag, seed=1, Nblock=3, nBoot=1000)
  sdmodel2.1 <- tsbootstrap(FUN = "ARywClassic",ts = tseries, order = lag, alpha1 = 1/0.54, alpha0 = 0, sigma_e = sigma_e,
                            seed = 1, Nblock=4, nBoot=1000)
  sdmodel2.2 <- tsbootstrap(FUN = "ARywClassic",ts = tseries, order = lag, alpha1 = 1/0.54, alpha0 = 0, sigma_e = sigma_e*2,
                            seed = 1, Nblock=4, nBoot=1000)
  sdmodel3.1 <- tsbootstrap(FUN  = "ARywMultiplicative",ts = tseries, order = lag, beta0 = 1/0.54, sigma_u = sigma_u,
                            seed = 1, Nblock=4, nBoot=1000)
  sdmodel3.2 <- tsbootstrap(FUN = "ARywMultiplicative",ts = tseries, order = lag, beta0 = 1/0.54, sigma_u = sigma_u*2,
                            seed = 1, Nblock=4, nBoot=1000)
  
  ### Print the Estimation Table  ####
  para <- unlist(list(model1,model2.1,model2.2,model3.1,model3.2))
  sd <- unlist(list(sdmodel1,sdmodel2.1,sdmodel2.2,sdmodel3.1,sdmodel3.2))
  
  Zvalue <- para/sd
  pvalue <- 2*(1-pt(abs(Zvalue),50))
  CIUpper <- para + qt(0.975, df = 50)*sd
  CILower <- para - qt(0.975, df = 50)*sd
  
  TableTemp <- cbind(para,sd)
  TableCI <- cbind(CILower,CIUpper)
  
  Table1qu_all <- data.frame(para= para, sd = sd, pvalue = pvalue)
  Table1qu_all <- round(Table1qu_all,3)
  Table1qu_all <-  data.frame(par = rep(c(rep("phi", lag+1),"sigmasqu"), 5 ), Table1qu_all)
  
  return(list(Table1 = Table1qu_all,
              model1 = model1,
              model2.1 = model2.1,
              model2.2 = model2.2,
              model3.1 = model3.1,
              model3.2 = model3.2,
              lag = lag, sigma_e = sigma_e, sigma_u = sigma_u))
}


diff_fit <- function (tseries, lag, sigma_e, sigma_u) {
  nodiff_fit(diff(tseries), lag, sigma_e, sigma_u)
 }

predictionInt_wrap <- function(lag, type, ...){
  if (type == "add"){
    if (lag == 1) { a <- predictionIntAR1_add(...)} else
      if (lag == 2) { a <- predictionIntAR2_add(...)} else
      { a <- predictionIntAR3_add(...)}
  }
  if (type == "mult"){
    if (lag == 1) { a <- predictionIntAR1_mult(...)} else
      if (lag == 2) { a <- predictionIntAR2_mult(...)} else
      { a <- predictionIntAR3_mult(...)}
  }
  if (type == "none"){
    if (lag == 1) { a <- predictionIntAR1(...)} else
      if (lag == 2) { a <- predictionIntAR2(...)} else
      { a <- predictionIntAR3(...)}
  }
  return(a)
}

prediction_nodiff <- function(tseries, Models, Province, fullseries, origindate = "2020-04-04"){
  laghere <- Models$lag
  sigma_e <- Models$sigma_e
  sigma_u <- Models$sigma_u
  
  ### Predicted the measurement error addressed time series
  
  mu <- mean(tseries)
  EstimateCase <- tseries[(length(tseries)-laghere):length(tseries)] * 0.54
  
  EstPred.Case1 <- prediction.AR( EstimateCase, 5, Models$model1, order = laghere)
  EstPred.Case1 <- EstPred.Case1[(laghere+1):length(EstPred.Case1)]
  EstPred.Case1UB <- EstPred.Case1 + predictionInt_wrap(lag = laghere, type = "none", 5, Models$model1)
  EstPred.Case1LB <- EstPred.Case1 - predictionInt_wrap(lag = laghere, type = "none", 5, Models$model1)
  Pe_h_1 <- predictionInt_wrap(lag = laghere, type = "none", 5, Models$model1)/1.96
  
  EstPred.Case2.1 <- prediction.AR( EstimateCase, 5, Models$model2.1, order = laghere)
  EstPred.Case2.1 <- EstPred.Case2.1[(laghere+1):length(EstPred.Case2.1)]
  EstPred.Case2.1UB <- EstPred.Case2.1 + predictionInt_wrap(lag = laghere, type = "add", 5, Models$model2.1, alpha1 = 1/0.54, sigma_e = sigma_e)
  EstPred.Case2.1LB <- EstPred.Case2.1 - predictionInt_wrap(lag = laghere, type = "add",5, Models$model2.1, alpha1 = 1/0.54, sigma_e = sigma_e)
  Pe_h_2.1 <- predictionInt_wrap(lag = laghere, type = "add", 5, Models$model2.1, alpha1 = 1/0.54, sigma_e = sigma_e)/1.96
  EstPred.Case2.1_Fitted <- FittedValue.AR(tseries * 0.54, 
                                           object = Models$model2.1, diff = F, order = laghere)
  
  
  EstPred.Case2.2 <- prediction.AR( EstimateCase, 5, Models$model2.2, order = laghere)
  EstPred.Case2.2 <- EstPred.Case2.2[(laghere+1):length(EstPred.Case2.2)]
  EstPred.Case2.2UB <- EstPred.Case2.2 + predictionInt_wrap(lag = laghere, type = "add", 5, Models$model2.2, alpha1 = 1/0.54, sigma_e = sigma_e *2)
  EstPred.Case2.2LB <- EstPred.Case2.2 - predictionInt_wrap(lag = laghere, type = "add", 5, Models$model2.2, alpha1 = 1/0.54, sigma_e = sigma_e *2)
  Pe_h_2.2 <- predictionInt_wrap(lag = laghere, type = "add", 5, Models$model2.2, alpha1 = 1/0.54, sigma_e = 0.1)/1.96
  EstPred.Case2.2_Fitted <- FittedValue.AR(tseries * 0.54, 
                                           object = Models$model2.2, diff = F, order = laghere)
  
  EstPred.Case1.1 <- prediction.AR( EstimateCase, 5, Models$model3.1, order = laghere)
  EstPred.Case1.1 <- EstPred.Case1.1[(laghere+1):length(EstPred.Case1.1)]
  EstPred.Case1.1UB <- EstPred.Case1.1 + predictionInt_wrap(lag = laghere, type = "mult", 5, Models$model3.1, mu = mu, sigma_u = sigma_u)
  EstPred.Case1.1LB <- EstPred.Case1.1 - predictionInt_wrap(lag = laghere, type = "mult", 5, Models$model3.1, mu = mu, sigma_u = sigma_u)
  Pe_h_3.1 <- predictionInt_wrap(lag = laghere, type = "mult", 5, Models$model3.1, mu = mu, sigma_u = sigma_u)/1.96
  EstPred.Case1.1_Fitted <- FittedValue.AR(tseries * 0.54, 
                                           object = Models$model3.1, diff = F, order = laghere)
  
  
  EstPred.Case1.2 <- prediction.AR( EstimateCase, 5, Models$model3.2, order = laghere)
  EstPred.Case1.2 <- EstPred.Case1.2[(laghere+1):length(EstPred.Case1.2)]
  EstPred.Case1.2UB <- EstPred.Case1.2 + predictionInt_wrap(lag = laghere, type = "mult", 5, Models$model3.2, mu = mu, sigma_u = sigma_u * 2)
  EstPred.Case1.2LB <- EstPred.Case1.2 - predictionInt_wrap(lag = laghere, type = "mult", 5, Models$model3.2, mu = mu, sigma_u = sigma_u * 2)
  Pe_h_3.2 <- predictionInt_wrap(lag = laghere, type = "mult", 5, Models$model3.2, mu = mu, sigma_u = sigma_u)/1.96
  EstPred.Case1.2_Fitted <- FittedValue.AR(tseries * 0.54, 
                                           object = Models$model3.2, diff = F, order = laghere)
  
  EstimateCasefull <- tseries * 0.54
  
  nn_model <- nnetar(tseries)
  nnetar_forcast <- forecast(EstimateCasefull, model = nn_model, h = 5, PI = TRUE, npaths = 100)
  nnetar_values <- c(EstimateCasefull[length(EstimateCasefull)], as.vector(nnetar_forcast$mean))
  
  ## rearrange the data into the format for plotting
  
  ### Plot ####
  
  # obtain the screeplot
  xlength <- length(EstPred.Case1)
  
  nfittedvalues <- length(EstPred.Case1.2_Fitted) + (laghere - 1)
  
  EstPredDataME <- data.frame(
    Day = rep(as.Date(1:xlength + nfittedvalues, origin = origindate),4),
    Estimate = c(EstPred.Case2.1,EstPred.Case2.2,EstPred.Case1.1,EstPred.Case1.2),
    EstimateUB = c(EstPred.Case2.1UB,EstPred.Case2.2UB,EstPred.Case1.1UB,EstPred.Case1.2UB),
    EstimateLB = c(EstPred.Case2.1LB,EstPred.Case2.2LB,EstPred.Case1.1LB,EstPred.Case1.2LB),
    MeasurementErrorDegree =rep(c("Mild","Moderate"),each=xlength,2) ,
    MeasurementErrorType = rep(c("Additive","Multiplicative"),each=xlength*2),
    MeasurementErrorTypeG = rep(c("Additive","Multiplicative"),each=xlength*2),
    TimeSeries = Province,
    Model = paste0("AR(",laghere,")"),
    AdjustReference = rep(fullseries[(length(fullseries)-5):length(fullseries)] * 0.54,4),
    Pe_h = c(Pe_h_2.1, Pe_h_2.2, Pe_h_3.1, Pe_h_3.2)
  )
  
  EstPredDataNaive <- data.frame(
    Day = as.Date(1:xlength + nfittedvalues, origin = origindate) ,
    Estimate = EstPred.Case1,
    EstimateUB = EstPred.Case1UB,
    EstimateLB = EstPred.Case1LB,
    MeasurementErrorTypeG =  "Naive",
    TimeSeries = Province,
    Pe_h = c(Pe_h_1),
    AdjustReference = fullseries[(length(fullseries)-5):length(fullseries)] * 0.54
  )
  
  EstPredDatannetar <- data.frame(
    Day = as.Date(1:xlength + nfittedvalues, origin = origindate) ,
    Estimate = nnetar_values,
    MeasurementErrorTypeG =  "NNAR",
    TimeSeries = Province
  )
  
  RealCase <- data.frame(
    Day = as.Date(1:length(fullseries), origin = origindate),
    Estimate =  fullseries,
    TimeSeries = Province,
    sgroup = "Reported Fatality"
  )
  
  
  AdjustedCase <- data.frame(
    Day = as.Date(1:length(fullseries), origin = origindate),
    Estimate =  fullseries * 0.54,
    TimeSeries = Province,
    sgroup = "Adjusted Fatality"
  )
  
  
  FittedValue <- data.frame(
    Day = rep(as.Date((laghere+1):(length(fullseries)-5), origin = origindate),4),
    Estimate =  c(EstPred.Case2.1_Fitted,
                  EstPred.Case2.2_Fitted,
                  EstPred.Case1.1_Fitted,
                  EstPred.Case1.2_Fitted),
    sgroup = "Fitted Fatality"
  )  
  
  return(list(Models = Models,
              Province = Province,
              EstPredDataME = EstPredDataME,
              EstPredDataNaive = EstPredDataNaive,
              EstPredDatannetar = EstPredDatannetar,
              RealCase = RealCase,
              AdjustedCase = AdjustedCase,
              FittedValue = FittedValue
  ))
}


prediction_diff <- function(tseries, Models, Province, fullseries, origindate = "2020-04-04"){
  laghere <- Models$lag
  sigma_e <- Models$sigma_e
  sigma_u <- Models$sigma_u
  
  ### Predicted the measurement error addressed time series
  
  mu <- mean(tseries)
  EstimateCase <- tseries[(length(tseries)-laghere):length(tseries)] * 0.54
  
  EstPred.Case1diff <- prediction.AR( diff(EstimateCase), 5, Models$model1, order = laghere)
  EstPred.Case1 <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case1diff[(laghere + 1):length(EstPred.Case1diff)]))
  EstPred.Case1UBdiff <- EstPred.Case1diff[laghere:length(EstPred.Case1diff)] + predictionInt_wrap(lag = laghere, type = "none", 5, Models$model1)
  EstPred.Case1UB <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case1UBdiff[2:length(EstPred.Case1UBdiff)])) 
  EstPred.Case1LBdiff <- EstPred.Case1diff[laghere:length(EstPred.Case1diff)] - predictionInt_wrap(lag = laghere, type = "none", 5, Models$model1)
  EstPred.Case1LB <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case1LBdiff[2:length(EstPred.Case1LBdiff)])) 
  Pe_h_1 <- predictionInt_wrap(lag = laghere, type = "none", 5, Models$model1)/1.96
  
  EstPred.Case2.1diff <- prediction.AR( diff(EstimateCase), 5, Models$model2.1, order = laghere)
  EstPred.Case2.1 <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case2.1diff[(laghere + 1):length(EstPred.Case2.1diff)]))
  EstPred.Case2.1UBdiff <- EstPred.Case2.1diff[laghere:length(EstPred.Case2.1diff)] + predictionInt_wrap(lag = laghere, type = "add", 5, Models$model2.1, alpha1 = 1/0.54, sigma_e = sigma_e)
  EstPred.Case2.1UB <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case2.1UBdiff[2:length(EstPred.Case2.1UBdiff)])) 
  EstPred.Case2.1LBdiff <- EstPred.Case2.1diff[laghere:length(EstPred.Case2.1diff)] - predictionInt_wrap(lag = laghere, type = "add", 5, Models$model2.1, alpha1 = 1/0.54, sigma_e = sigma_e)
  EstPred.Case2.1LB <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case2.1LBdiff[2:length(EstPred.Case2.1LBdiff)])) 
  Pe_h_2.1 <- predictionInt_wrap(lag = laghere, type = "add", 5, Models$model2.1, alpha1 = 1/0.54, sigma_e = sigma_e)/1.96
  EstPred.Case2.1_Fitted <- FittedValue.AR(fullseries[1:(length(fullseries)-5)] * 0.54, 
                                           object = Models$model2.1, diff = T, order = laghere)
  
  EstPred.Case2.2diff <- prediction.AR( diff(EstimateCase), 5, Models$model2.2, order = laghere)
  EstPred.Case2.2 <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case2.2diff[(laghere + 1):length(EstPred.Case2.2diff)]))
  EstPred.Case2.2UBdiff <- EstPred.Case2.2diff[laghere:length(EstPred.Case2.2diff)] + predictionInt_wrap(lag = laghere, type = "add", 5, Models$model2.2, alpha1 = 1/0.54, sigma_e = sigma_e *2)
  EstPred.Case2.2UB <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case2.2UBdiff[2:length(EstPred.Case2.2UBdiff)])) 
  EstPred.Case2.2LBdiff <- EstPred.Case2.2diff[laghere:length(EstPred.Case2.2diff)] - predictionInt_wrap(lag = laghere, type = "add", 5, Models$model2.2, alpha1 = 1/0.54, sigma_e = sigma_e *2)
  EstPred.Case2.2LB <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case2.2LBdiff[2:length(EstPred.Case2.2LBdiff)])) 
  Pe_h_2.2 <- predictionInt_wrap(lag = laghere, type = "add", 5, Models$model2.2, alpha1 = 1/0.54, sigma_e = sigma_e *2)/1.96
  EstPred.Case2.2_Fitted <- FittedValue.AR(fullseries[1:(length(fullseries)-5)] * 0.54, 
                                           object = Models$model2.2, diff = T, order = laghere)
  
  EstPred.Case3.1diff <- prediction.AR( diff(EstimateCase), 5, Models$model3.1, order = laghere)
  EstPred.Case3.1 <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case3.1diff[(laghere + 1):length(EstPred.Case3.1diff)]))
  EstPred.Case3.1UBdiff <- EstPred.Case3.1diff[laghere:length(EstPred.Case3.1diff)] + predictionInt_wrap(lag = laghere, type = "mult", 5, Models$model3.1, mu = mu, sigma_u = sigma_u)
  EstPred.Case3.1UB <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case3.1UBdiff[2:length(EstPred.Case3.1UBdiff)])) 
  EstPred.Case3.1LBdiff <- EstPred.Case3.1diff[laghere:length(EstPred.Case3.1diff)] - predictionInt_wrap(lag = laghere, type = "mult", 5, Models$model3.1, mu = mu, sigma_u = sigma_u)
  EstPred.Case3.1LB <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case3.1LBdiff[2:length(EstPred.Case3.1LBdiff)])) 
  Pe_h_3.1 <- predictionInt_wrap(lag = laghere, type = "mult", 5, Models$model3.1, mu = mu, sigma_u = sigma_u)/1.96
  EstPred.Case3.1_Fitted <- FittedValue.AR(fullseries[1:(length(fullseries)-5)] * 0.54, 
                                           object = Models$model3.1, diff = T, order = laghere)
  
  EstPred.Case3.2diff <- prediction.AR( diff(EstimateCase), 5, Models$model3.2, order = laghere)
  EstPred.Case3.2 <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case3.2diff[(laghere + 1):length(EstPred.Case3.2diff)]))
  EstPred.Case3.2UBdiff <- EstPred.Case3.2diff[laghere:length(EstPred.Case3.2diff)] + predictionInt_wrap(lag = laghere, type = "mult", 5, Models$model3.2, mu = mu, sigma_u = sigma_u * 2)
  EstPred.Case3.2UB <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case3.2UBdiff[2:length(EstPred.Case3.2UBdiff)])) 
  EstPred.Case3.2LBdiff <- EstPred.Case3.2diff[laghere:length(EstPred.Case3.2diff)] - predictionInt_wrap(lag = laghere, type = "mult", 5, Models$model3.2, mu = mu, sigma_u = sigma_u * 2)
  EstPred.Case3.2LB <- cumsum(c(EstimateCase[length(EstimateCase)],EstPred.Case3.2LBdiff[2:length(EstPred.Case3.2LBdiff)])) 
  Pe_h_3.2 <- predictionInt_wrap(lag = laghere, type = "mult", 5, Models$model3.2, mu = mu, sigma_u = sigma_u * 2)/1.96
  EstPred.Case3.2_Fitted <- FittedValue.AR(fullseries[1:(length(fullseries)-5)] * 0.54, 
                                           object = Models$model3.2, diff = T, order = laghere)
  
  
  
  EstimateCasefull <- tseries * 0.54
  
  nn_model <- nnetar(tseries)
  nnetar_forcast <- forecast(EstimateCasefull, model = nn_model, h = 5, PI = TRUE, npaths = 100)
  nnetar_values <- c(EstimateCasefull[length(EstimateCasefull)], as.vector(nnetar_forcast$mean))
  
  ## rearrange the data into the format for plotting
  
  ### Plot ####
  
  # obtain the screeplot
  xlength <- length(EstPred.Case1)
  
  nfittedvalues <- length(EstPred.Case3.2_Fitted) + laghere 
  
  EstPredDataME <- data.frame(
    Day = rep(as.Date(1:xlength + nfittedvalues, origin = origindate),4),
    Estimate = c(EstPred.Case2.1,EstPred.Case2.2,EstPred.Case3.1,EstPred.Case3.2),
    EstimateUB = c(EstPred.Case2.1UB,EstPred.Case2.2UB,EstPred.Case3.1UB,EstPred.Case3.2UB),
    EstimateLB = c(EstPred.Case2.1LB,EstPred.Case2.2LB,EstPred.Case3.1LB,EstPred.Case3.2LB),
    MeasurementErrorDegree =rep(c("Mild","Moderate"),each=xlength,2) ,
    MeasurementErrorType = rep(c("Additive","Multiplicative"),each=xlength*2),
    MeasurementErrorTypeG = rep(c("Additive","Multiplicative"),each=xlength*2),
    TimeSeries = Province,
    Model = paste0("AR(",laghere,")"),
    AdjustReference = rep(fullseries[(length(fullseries)-5):length(fullseries)] * 0.54,4),
    Pe_h = c(Pe_h_2.1, Pe_h_2.2, Pe_h_3.1, Pe_h_3.2)
  )
  
  EstPredDataNaive <- data.frame(
    Day = as.Date(1:xlength + nfittedvalues, origin = origindate) ,
    Estimate = EstPred.Case1,
    EstimateUB = EstPred.Case1UB,
    EstimateLB = EstPred.Case1LB,
    MeasurementErrorTypeG =  "Naive",
    TimeSeries = Province,
    Pe_h = c(Pe_h_1),
    AdjustReference = fullseries[(length(fullseries)-5):length(fullseries)] * 0.54
  )
  
  EstPredDatannetar <- data.frame(
    Day = as.Date(1:xlength + nfittedvalues, origin = origindate) ,
    Estimate = nnetar_values,
    MeasurementErrorTypeG =  "NNAR",
    TimeSeries = Province
  )
  
  RealCase <- data.frame(
    Day = as.Date(1:length(fullseries), origin = origindate),
    Estimate =  fullseries,
    TimeSeries = Province,
    sgroup = "Reported Fatality"
  )
  
  
  AdjustedCase <- data.frame(
    Day = as.Date(1:length(fullseries), origin = origindate),
    Estimate =  fullseries * 0.54,
    TimeSeries = Province,
    sgroup = "Adjusted Fatality"
  )
  
  
  FittedValue <- data.frame(
    Day = rep(as.Date((laghere+2):(length(fullseries)-5), origin = origindate),4),
    Estimate =  c(EstPred.Case2.1_Fitted,
                  EstPred.Case2.2_Fitted,
                  EstPred.Case3.1_Fitted,
                  EstPred.Case3.2_Fitted),
    sgroup = "Fitted Fatality"
  )  
  
  return(list(Models = Models,
              Province = Province,
              EstPredDataME = EstPredDataME,
              EstPredDataNaive = EstPredDataNaive,
              EstPredDatannetar = EstPredDatannetar,
              RealCase = RealCase,
              AdjustedCase = AdjustedCase,
              FittedValue = FittedValue
  ))
}




plot_results <- function (PredictionResults,defnumber, yrange = c(0,3), diff = F){
  wdiff <- ifelse(diff, "diff", "nodiff")
  pdf(file = paste0("output/revision/2021-01-15-",PredictionResults$Province,"-def",
                    defnumber,"-AR(",PredictionResults$Models$lag,")-", wdiff,".pdf"), height = 10 , width = 18)
  strip_color <- "#0A1D37"
  p1 <- ggplot(PredictionResults$EstPredDataME)  +
    geom_line(aes(x=Day, y=Estimate, alpha = sgroup), data = PredictionResults$RealCase, col="#150718",size=1)+
    geom_line(aes(x=Day, y=Estimate, alpha = sgroup), data = PredictionResults$AdjustedCase, col="#6D882B",size=1)+
    geom_line(aes(x=Day, y=Estimate, alpha = sgroup), data = PredictionResults$FittedValue, col="#5D6D7E",size=1, linetype = 3)+
    geom_ribbon(aes(x=Day,ymin = EstimateLB, ymax = EstimateUB, fill = MeasurementErrorTypeG), alpha=0.12) +
    geom_ribbon(aes(x=Day,ymin = EstimateLB, ymax = EstimateUB, fill = MeasurementErrorTypeG), alpha=0.12) +
    geom_line(aes(x=Day, y=Estimate, color = MeasurementErrorTypeG),size=1.2,data = PredictionResults$EstPredDataNaive, alpha= 0.7) +
    geom_point(aes(x=Day, y=Estimate, color = MeasurementErrorTypeG, shape= MeasurementErrorTypeG),size=3,data = PredictionResults$EstPredDataNaive) +
    geom_line(aes(x=Day, y=Estimate, color = MeasurementErrorTypeG),size=1.2,data = PredictionResults$EstPredDatannetar, alpha= 0.7) +
    geom_point(aes(x=Day, y=Estimate, color = MeasurementErrorTypeG, shape= MeasurementErrorTypeG),size=3,data = PredictionResults$EstPredDatannetar) +
    geom_line(aes(x=Day, y=Estimate, color = MeasurementErrorTypeG, linetype=MeasurementErrorDegree), size=1.2)+
    geom_point(aes(x=Day, y=Estimate, color = MeasurementErrorTypeG, shape= MeasurementErrorTypeG),size=2) + 
    # scale_color_manual(values = c("#469BEC", "#C4878C","#AC9344"))  +
    scale_color_manual(values = c("#516888", "#802729","#F3AE6D","#8CBEB1"))  +
    scale_fill_manual(values = c("#516888", "#802729","#F1E3B6"))  +  
    scale_shape_manual(values = c(19, 19, 10, 17))  +
    labs(x = "Day", y = "Fatality Rate (%)") +
    scale_alpha_manual(values=c(0.7, 1, 0.7),
                       guide=guide_legend(override.aes = list(colour=c("#6D882B","#5D6D7E","#150718"), linetype = c(1,3,1), alpha = c(1,1,1)))) +
    scale_linetype_manual(values=c(3, 5),guide = FALSE) +
    labs(color = "Methods", alpha = "Reference Time Series",  shape = "Methods", fill = "95% Prediction Interval") +
    geom_vline(xintercept=c(as.Date("2021-01-10")), linetype=2)+
    coord_cartesian(xlim = c(as.Date("2020-12-10"), NA), ylim = yrange) +
    facet_grid(MeasurementErrorDegree ~ MeasurementErrorType) +
    theme_bw() +
    theme(strip.background =element_rect(fill=strip_color,color=strip_color))+ # #535b44
    theme(strip.text = element_text(colour = 'white')) +
    theme(panel.border = element_rect(colour = strip_color))+
    theme(text=element_text(size=23, family="URWHelvetica"))
  
  print(p1)
  
  
  dev.off()
  
}


Evaluations_data <- function(PredictionResults){
  Pe_h_est_pro <- PredictionResults$EstPredDataME %>% 
    mutate(Pe_h_est = (AdjustReference-Estimate)^2) %>% 
    mutate(Pe_h = Pe_h^2) %>% 
    select("Day","Pe_h_est", "Pe_h", "MeasurementErrorDegree","MeasurementErrorTypeG")
  
  Pe_h_est_naive <- PredictionResults$EstPredDataNaive %>% 
    mutate(Pe_h_est = (AdjustReference-Estimate)^2) %>% 
    mutate(Pe_h = Pe_h^2) %>% 
    mutate(MeasurementErrorDegree = NA) %>%
    select("Day","Pe_h_est", "Pe_h", "MeasurementErrorDegree","MeasurementErrorTypeG")
  
  Pe_h_est <- rbind(Pe_h_est_naive, Pe_h_est_pro)
  
  Pe_h_est_wide <- Pe_h_est %>% select("Day","Pe_h_est","MeasurementErrorTypeG", "MeasurementErrorDegree") %>%
    spread(Day, c("Pe_h_est")) %>% mutate(sums = rowSums(.[4:8]))  %>% mutate(Model = paste0("AR(", PredictionResults$Models$lag,")"))
  
  Pe_h_wide <- Pe_h_est %>% select("Day","Pe_h","MeasurementErrorTypeG", "MeasurementErrorDegree") %>%
    spread(Day, c("Pe_h")) %>% mutate(sums = rowSums(.[4:8]))  %>% mutate(Model = paste0("AR(", PredictionResults$Models$lag,")"))
  return(list(Pe_h_est_pro = Pe_h_est_pro,
              Pe_h_est_naive = Pe_h_est_naive,
              Pe_h_est = Pe_h_est,
              Pe_h_est_wide = Pe_h_est_wide,
              Pe_h_wide = Pe_h_wide))
}


find_optimal_lag <- function(ts){
  AIC_current <- Inf
  CONTINUE <- T
  initi <- 1
  while (CONTINUE) {
    fit <- Arima(ts, order= c(initi, 0, 0))
    cat(paste0("Lag ", initi, " AIC:", fit$aic))
    if (fit$aic > AIC_current) {CONTINUE <- F}
    else {
      initi <- initi + 1
      AIC_current <- fit$aic }
  }
  return(initi - 1)
}


get_Peh_table_prov <- function(Pe_h_est_list){
  Table <- NULL
  for (i in 1:length(Pe_h_est_list)) {
    Table <- rbind(Table, cbind(Pe_h_est_list[[i]]$Pe_h_est_wide[c(1,2,10,4:9)],Pe_h_est_list[[i]]$Pe_h_wide[c(4:9)]))}
  
  colnames(Table) <- c("MeasurementErrorTypeG", "MeasurementErrorDegree", "Model", "EP2021-01-10", "EP2021-01-11",
                       "EP2021-01-12", "EP2021-01-13", "EP2021-01-14", "sums1", "2021-01-10",
                       "2021-01-11", "2021-01-12", "2021-01-13", "2021-01-14", "sums2")
  
  Table <- Table %>% 
    mutate_if(is.numeric, round, digits = 3)
  
  return(Table)
}


produce_pacf_plot <- function(tseries, Province, defnumber){
  pdf(file = paste0("output/revision/PACF-",Province,"-def", defnumber,".pdf"), height = 6 , width = 8)
  pacfplot <- pacf(tseries, main=paste0("PACF plot: ", Province," (Definition ", defnumber, ")"), cex.main = 5)
  print(pacfplot)
  dev.off()
  
  pdf(file = paste0("output/revision/ACF-",Province,"-def", defnumber,".pdf"), height = 6 , width = 8)
  pacfplot <- acf(tseries, main=paste0("ACF plot: ", Province," (Definition ", defnumber, ")"), cex.main = 5)
  print(pacfplot)
  dev.off()
}
