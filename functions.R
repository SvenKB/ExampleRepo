#Circular mediation
# Load necessary packages
library(circular)
library(circglmbayes)
library(coda)
library(boot)
library(parallel)



# Difference in coefficients

CircMed_Diff <- function(dt,ind, predictor = "x", mediator = "m", outcome = "y") {
  
  # Dataset of this iteration
  dat <- dt[ind,]
  
  # Standardize predictors
  x <- dat[,predictor]
  m <- dat[,mediator]
  y <- dat[,outcome]
  
  x <- (x-mean(x))/sqrt(var(x))
  m <- (m-mean(m))/sqrt(var(m))
  
  # Create predictor matrix
  predictors <- cbind(x,m)
  # Prepare outcome 
  outcome <- circular:::as.circular(y)
  
  boot_dt <- data.frame(predictors,y=outcome)
  # Models
  mediator_model <- lm(m~x, data=dat)
  mediated_model <- circular:::lm.circular.cl(y = outcome ,x = predictors , init = c(0,0))
  direct_model <- circular:::lm.circular.cl(y = outcome,x = predictors[,1], init = 0)
  
  # Coefficients
  a <- mediator_model$coefficients[[2]]
  b <- mediated_model$coefficients[[2]]
  c <- mediated_model$coefficients[[1]]
  c_tilde <- direct_model$coefficients
  
  # Calculate effects
  direct_effect <- c
  total_effect <- c_tilde
  indirect_effect <- c_tilde - c

  # Prepare output
  kappa <- mediated_model$kappa
  output <- c(total_effect,direct_effect,indirect_effect,a,b,kappa)
  names(output) <- c("Total","Direct","Indirect","a","b","Residual Kappa")
  
  return(output)
}


CircMed_Product <- function(dt,ind, predictor = "x", mediator = "m", outcome = "y") {
  
  # Dataset of this iteration
  dat <- dt[ind,]
  
  # Standardize predictors
  x <- dat[,predictor]
  m <- dat[,mediator]
  y <- dat[,outcome]
  
  x <- (x-mean(x))/sqrt(var(x))
  m <- (m-mean(m))/sqrt(var(m))
  
  # Create predictor matrix
  predictors <- cbind(x,m)
  # Prepare outcome 
  outcome <- as.circular(y)
  
  # Models
  dat <- data.frame(x,m)
  mediator_model <- lm(m~x, data=dat)
  mediated_model <- circular:::lm.circular.cl(y = outcome ,x = predictors , init = c(0,0))
  total_model <- circular:::lm.circular.cl(y = outcome,x = predictors[,1], init = 0)
  
  # Coefficients
  a <- mediator_model$coefficients[[2]]
  b <- mediated_model$coefficients[[2]]
  c <- mediated_model$coefficients[[1]]
  c_tilde <- total_model$coefficients
  
  
  # Calculate effects
  indirect_effect <- a*b
  direct_effect <- c
  total_effect <- c_tilde
  
  # Prepare output
  kappa <- mediated_model$kappa
  output <- c(total_effect,direct_effect,indirect_effect,a,b,kappa)
  names(output) <- c("Total","Direct","Indirect","a","b","Residual kappa")
  
  return(output)
  
}




CircMed_Reparameter <- function(dt,ind, predictor = "x", mediator = "m", outcome = "y") {
  
  # Dataset of this iteration
  dat <- dt[ind,]
  
  # Standardize predictors
  x <- dat[,predictor]
  m <- dat[,mediator]
  y <- dat[,outcome]
  
  x <- (x-mean(x))/sqrt(var(x))
  m <- (m-mean(m))/sqrt(var(m))
  
  # x-residualize M
  dat <- data.frame(x,m)
  resid_model <- lm(m~x, data = dat)
  m_tilde <- resid_model$residuals
  
  # Create predictor matrix
  predictors <- cbind(x,m_tilde)
  # Prepare outcome 
  outcome <- as.circular(y)
  
  # Models
  # Direct model
  mediated_model <- circular:::lm.circular.cl(y = outcome ,x = predictors , init = c(0,0))
  # a-path model
  mediator_model <- lm(m~x,data=dat)
  
  #Coefficients
  a <- mediator_model$coefficients[[2]]
  b <- mediated_model$coefficients[[2]]
  # Calculate effects
  total_effect <- mediated_model$coefficients[[1]]
  indirect_effect <- mediator_model$coefficients[[2]]*mediated_model$coefficients[[2]]
  direct_effect <- total_effect-indirect_effect
  
  # Prepare output
  kappa <- mediated_model$kappa
  output <- c(total_effect,direct_effect,indirect_effect,a,b,kappa)
  names(output) <- c("Total","Direct","Indirect","a","b","Resid. Kappa")
  return(output)
}




CircMed_Bayes_Diff <- function(dat) {
 
  # Standardize predictors
  x <- dat[,1]
  m <- dat[,2]
  y <- dat[,3]
  
  x <- (x-mean(x))/sqrt(var(x))
  m <- (m-mean(m))/sqrt(var(m))
  
  y <- as.circular(y)
  
  # Create dataframe
  data <- data.frame(x,m,y)
  
  
  # Models
  pred_med <- lm(m~x, data=data)
  a <- rnorm(50000,pred_med$coefficients[2], summary(pred_med)$coefficients[2,2])
  total_model    <- circGLM(y ~ x, data = data, Q = 50000)
  mediated_model  <- circGLM(y ~ x+m, data = data, Q = 50000)
  
  # The posterior sample of mediation effects using difference method. 
  mediation_sample_difference <- cbind(total    = total_model$bt_chain[,1], 
                                       direct   = mediated_model$bt_chain[,1],
                                       indirect = total_model$bt_chain[,1] - mediated_model$bt_chain[,1],
                                       b =  mediated_model$bt_chain[,2],
                                       a = a)
  
  # Obtain a summary of effects
  mcmcsum <- summary(mcmc(mediation_sample_difference)) 
  
  # Combine into a table. 
  #list(cbind(mcmcsum$statistics, mcmcsum$quantiles),kappa=mediated_model$kp_mean)
  #output <- list(mcmcsum$statistics[1,1],
  #               mcmcsum$statistics[2,1],
  #               mcmcsum$statistics[3,1],
  #               mediated_model$kp_mean)
  #names(output) <- c("Total","Direct","Indirect","Resid. Kappa")
  output <- list(mcmcsum$statistics, mcmcsum$quantiles,kappa=mediated_model$kp_mean)
  
  return(output)
}

CircMed_Bayes_Product <- function (dat) {
  
  # Standardize predictors
  x <- dat[,1]
  m <- dat[,2]
  y <- dat[,3]
  
  x <- (x-mean(x))/sqrt(var(x))
  m <- (m-mean(m))/sqrt(var(m))
  
  y <- as.circular(y)
  
  # Create dataframe
  data <- data.frame(x,m,y)
  
  ## Models
  
  # Predictor-Mediator model
  #pred_med <- brm(m~x, data=data, warmup = 1000, iter = 51000, chains = 1)
  pred_med <- lm(m~x, data=data)
  a <- rnorm(50000,pred_med$coefficients[2], summary(pred_med)$coefficients[2,2])
  # Mediated model
  mediated_model  <- circGLM(y ~ x+m, data = data, Q = 50000)
  
  # Total model
  total_model <- circGLM(y ~ x, data = data, Q = 50000)
  
  mediation_sample_product <- cbind(total    = total_model$bt_chain[,1], 
                                    direct   = mediated_model$bt_chain[,1],
                                    indirect = mediated_model$bt_chain[,2]*a,
                                    a = a,
                                    b = mediated_model$bt_chain[,2])
  
  # Obtain a summary of effects
  mcmcsum <- summary(mcmc(mediation_sample_product)) 
  
  # Combine into a table. 
  #list(cbind(mcmcsum$statistics, mcmcsum$quantiles),kappa=mediated_model$kp_mean)
  #output <- list(mcmcsum$statistics[1,1],
  #               mcmcsum$statistics[2,1],
  #               mcmcsum$statistics[3,1],
  #               mediated_model$kp_mean)
  #names(output) <- c("Total","Direct","Indirect","Resid. Kappa")
  output <- list(mcmcsum$statistics, mcmcsum$quantiles,kappa=mediated_model$kp_mean)
  return(output)
}


simData <- function(a,b,c,n) {
 
  linkfun   = function(x) 2 * atan(x)

  x <- rnorm(n,0,1)
  m <- rnorm(n,(a*x),1)
  y <- rep(0,n)

  beta <- c(c,b)
  pred <- cbind(x,m)
 
  con <- linkfun(apply(pred, 1, "%*%", beta))
  
  y_pred <- 1+con 
  err <- rvmc(n,0,5)
  y <- y_pred + err
  y <- as.circular(y)
  
  data <- data.frame(x,m,y)
  
  return(data)
}

mediationBootstrap <- function(dt, fun, R = 100, probs = c(.025,.25,.5,.75,.975)) {
  
  
  circmedboot <- boot(data=dt,fun, R = R, stype = "i")
  
  
  ps <- apply(circmedboot$t, 2L, function(x) mean(x <= 0))
  
  
  res <- list(tab = cbind(original     = circmedboot$t0,
                          bias         = apply(circmedboot$t, 2L, mean, na.rm = TRUE) - circmedboot$t0,
                          "std. error" = sqrt(apply(circmedboot$t, 2L, function(t.st) var(t.st[!is.na(t.st)]))),
                          t(apply(circmedboot$t, 2L, function(x) quantile(x, probs = probs))),
                          "one-sided p-value"    = ifelse(circmedboot$t0 > 0, ps, 1 - ps),
                          "two-sided p-value"    = ifelse(circmedboot$t0 > 0, ps*2, (1 - ps) * 2 )),
              bootsam = circmedboot$t,
              bootobj = circmedboot)
  res
}




#### Circular-Linear
## ------------------------------------------------------------------------
# Calculate the standardized path coefficient between x and y given z.
path_coef <- function(y, x, z) unname(coef(lm(scale(y) ~ scale(x) + z))[2])

## ------------------------------------------------------------------------
# Euclidean norm for later use. 
l2norm <- function(x) sqrt(x[1]^2 + x[2]^2)

# Compute the circular-linear correlation model, with linear predictor x,
# circular outcome th and linear mediator z.
clcor_mediation <- function(th,x,z) {
  
  # Rotation by second eigenvector.
  evec2 <- eigen(cov(cbind(cos(th),  sin(th))))$vectors[, 2]
  rot   <- atan2(evec2[2], evec2[1]) 
  
  # Vector of angles with uncorrelated sine and cosine. 
  th_rot <- th - rot
  
  # First get the results as (cosine, sine) vectors. 
  # Total effect.
  tv <- c(cor(cos(th_rot), x), cor(sin(th_rot), x))
  
  # Direct effect
  dv <- c(path_coef(cos(th_rot), x, z), path_coef(sin(th_rot), x, z))
  
  # Indirect effect.
  a <- cor(x, z)
  p_cz <- path_coef(cos(th_rot), z, x)
  p_sz <- path_coef(sin(th_rot), z, x)
  iv <- c(a * p_cz, a * p_sz)
  
  list(effects = c(total    = l2norm(tv), 
                   direct   = l2norm(dv), 
                   indirect = l2norm(iv)), 
       vectors = cbind(total    = tv,
                       direct   = dv,
                       indirect = iv))
}


computeClcorMediation <- function(data, inds, outcome="y", predictor="x", mediators="m") {
  listres <- clcor_mediation(th = data[inds, outcome], 
                             z  = data[inds, mediators], 
                             x  = data[inds, predictor])
  res        <- unlist(listres)
  tdi        <- c("Total", "Direct", "Indirect")
  names(res) <- c(tdi, paste0(rep(tdi, each = 2), c("_cos", "_sin")))
  res
}


clcorMediationBootstrap <- function(data, outcome="y", predictor="x", mediators="m", 
                                    R = 1000, probs = c(.025,.25,.5,.75,.975), ...) {
  
  
  # Compute the bootstrap using the boot package.
  suppressWarnings({
    clcormedboot <- boot(data, computeClcorMediation, R = R, 
                         outcome   = outcome, 
                         predictor = predictor, 
                         mediators = mediators, ...)
  })
  
  # Obtain the p-values for t0 > 0, which will be taken as 1 - p for t0 < 0.
  ps <- apply(clcormedboot$t, 2L, function(x) mean(x <= 0))
  
  # Obtain result
  res <- list(tab = cbind(original     = clcormedboot$t0,
                          bias         = apply(clcormedboot$t, 2L, mean, na.rm = TRUE) - clcormedboot$t0,
                          "std. error" = sqrt(apply(clcormedboot$t, 2L, function(t.st) var(t.st[!is.na(t.st)]))),
                          t(apply(clcormedboot$t, 2L, function(x) quantile(x, probs = probs))),
                          "one-sided p-value"    = ifelse(clcormedboot$t0 > 0, ps, 1 - ps),
                          "two-sided p-value"    = ifelse(clcormedboot$t0 > 0, ps*2, (1 - ps) * 2 )),
              bootsam = clcormedboot$t,
              bootobj = clcormedboot)
  
  class(res) <- c("clcorMedBoot", class(res))
  res
}

## ---------------------------------------------------------------

## RDS
loadDatasets <- function(truen,truea,trueb,truec,nsim) {
  
  #Prepare all possible designs
  Alldesigns <- expand.grid(a=truea,b=trueb,c=truec,n=truen,stringsAsFactors=FALSE)
  nonexistendDesigns <- 0
  data <- list()
  for (i in 1:nrow(Alldesigns)) {
    design <- Alldesigns[i,]
    
    curr_a <- design[,1]
    curr_b <- design[,2]
    curr_c <- design[,3]
    curr_n <- design[,4]
    
    # Prepare loading datasets
    DirName <- paste0(getwd(),
                      "/Data/Datasets_",
                      "n=", curr_n,
                      "a=", curr_a,
                      "b=", curr_b,
                      "c=", curr_c)
    DesignName <- paste0("n=", curr_n,
                         "a=", curr_a,
                         "b=", curr_b,
                         "c=", curr_c)
    dat <- list()
    # Load all datasets per design
    for (j in 1:nsim) {
      
      filename <- paste0(DirName,"/nr",j,".RDS")
      filenumber <- paste0("nr",j)
      
      if (file.exists(filename)) {
        dat[[filenumber]] <- readRDS(filename)
        names(dat[[filenumber]]) <- c("x","m","y")
        
      } else {
        nonexistendDesigns <- nonexistendDesigns + 1
      }
      
    }
    # Store list of datasets per design in a list
    data[[DesignName]] <- dat
  }
  if (nonexistendDesigns > 0) {
    cat("\n[Data loading: ", nonexistendDesigns, "/", nsim*nrow(Alldesigns),
        " datasets did not exist.]\n")
  }
  return(data)
}
## CSV
loadDatasets <- function(truen,truea,trueb,truec,nsim) {
  
  #Prepare all possible designs
  Alldesigns <- expand.grid(a=truea,b=trueb,c=truec,n=truen,stringsAsFactors=FALSE)
  nonexistendDesigns <- 0
  data <- list()
  for (i in 1:nrow(Alldesigns)) {
    design <- Alldesigns[i,]
    
    curr_a <- design[,1]
    curr_b <- design[,2]
    curr_c <- design[,3]
    curr_n <- design[,4]
    
    # Prepare loading datasets
    DirName <- paste0(getwd(),
                      "/Data/Datasets_",
                      "n=", curr_n,
                      "a=", curr_a,
                      "b=", curr_b,
                      "c=", curr_c)
    DesignName <- paste0("n=", curr_n,
                         "a=", curr_a,
                         "b=", curr_b,
                         "c=", curr_c)
    dat <- list()
    # Load all datasets per design
    for (j in 1:nsim) {
      
      filename <- paste0(DirName,"/nr",j,".csv")
      filenumber <- paste0("nr",j)
      
      if (file.exists(filename)) {
        dat[[filenumber]] <- read.csv(filename,header = FALSE)
        names(dat[[filenumber]]) <- c("x","m","y")
        
      } else {
        nonexistendDesigns <- nonexistendDesigns + 1
      }
      
    }
    # Store list of datasets per design in a list
    data[[DesignName]] <- dat
  }
  if (nonexistendDesigns > 0) {
    cat("\n[Data loading: ", nonexistendDesigns, "/", nsim*nrow(Alldesigns),
        " datasets did not exist.]\n")
  }
  return(data)
}


saveDatasets <- function(truen,truea,trueb,truec,nsim, seed = 140689) {
  set.seed(seed)
  
  # prepare all possible designs
  Alldesigns <- expand.grid(a=truea,b=trueb,c=truec,n=truen,stringsAsFactors=FALSE)
  existingDesigns <- 0
  
  for (i in 1:nrow(Alldesigns)) {
    design <- Alldesigns[i,]
    
    curr_a <- design[,1]
    curr_b <- design[,2]
    curr_c <- design[,3]
    curr_n <- design[,4]
    
    # Prepare saving datasets
    DirName <- paste0(getwd(),
                      "/Data/Datasets_",
                      "n=", curr_n,
                      "a=", curr_a,
                      "b=", curr_b,
                      "c=", curr_c)
    dir.create(DirName, showWarnings = FALSE)
    
    # Simulate n datasets per design
    for (j in 1:nsim) {
      
      filename <- paste0(DirName,"/nr",j,".csv")
      
      if (!file.exists(filename)) {
        dat <- simData(a=curr_a,b=curr_b,c=curr_c,n=curr_n)
        write.table(dat, filename, sep = ",", row.names=FALSE, col.names=FALSE)
      } else {
        existingDesigns <- existingDesigns + 1
      }
      
    }
    
  }
  if (existingDesigns > 0) {
    cat("\n[Data generation: ", existingDesigns, "/", nsim*nrow(Alldesigns),
        " datasets already existed.]\n")
  }
}

###################################################################################################
###################################################################################################
###################################################################################################