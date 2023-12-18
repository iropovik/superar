# Load libraries (and install if not installed already)
list.of.packages <- c("car", "reshape", "tidyverse", "tools", "tidyr", "psych", "metafor", "meta", "dmetar", "esc", "lme4", "ggplot2", "knitr", "puniform", "kableExtra", "lmerTest", "pwr", "Amelia", "multcomp", "magrittr", "weightr", "clubSandwich", "ddpcr", "poibin", "compute.es", "readxl", "robvis", "RoBMA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load required libraries
lapply(list.of.packages, require, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)

# Determine the alpha level / use of one-tailed vs two-tailed test for p-uniform and PET-PEESE
if (test == "one-tailed") {
  alpha <- .10
} else if (test == "two-tailed") {
  alpha <- .05
}

# Meta-analysis -----------------------------------------------------------

# Meta analysis run on a filtered dataset

# Custom robust (RVE) multivariate RE meta-analytic model (using the CHE working model)
# Needs specific naming of ES, variances and data on clustering; yi = yi, vi = vi, paperID, effectID
# Using n/(n-p) small-sample correction for RVE SEs
rmaCustom <- function(data = NA){
  viMatrix <- data %$% impute_covariance_matrix(vi, cluster = paperID, r = rho)
  rmaObjectModBasedSE <- rma.mv(yi = yi, V = viMatrix, data = data, method = rmaMethod, random = ~ 1|paperID/effectID, sparse = TRUE)
  rmaObject <- robust.rma.mv(rmaObjectModBasedSE, cluster = data$paperID)
  return(list("RMA.MV object with RVE SEs with n/(n-p) small-sample correction" = rmaObject, 
              "RMA.MV object without cluster robust SEs" = rmaObjectModBasedSE))
}

# 95% prediction interval -------------------------------------------------
pi95 <- function(rmaObject = NA){
  pi95Out <- c("95% PI LB" = round(predict.rma(rmaObject[[1]])$cr.lb, 3), "95% PI UB" = round(predict.rma(rmaObject[[1]])$cr.ub, 3))
  pi95Out
}

# Heterogeneity -----------------------------------------------------------

heterogeneity <- function(rmaObject = NA){
  
  # Total heterogeneity - tau
  tau <- sqrt(sum(rmaObject[[1]]$sigma2))
  
  # I^2
  W <- diag(1/rmaObject[[1]]$vi)
  X <- model.matrix(rmaObject[[1]])
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2<- 100 * sum(rmaObject[[1]]$sigma2) / (sum(rmaObject[[1]]$sigma2) + (rmaObject[[1]]$k - rmaObject[[1]]$p)/sum(diag(P)))
  
  # Separate estimates of between- and within-cluster heterogeneity
  BW.hetero <- round(100 * rmaObject[[1]]$sigma2 / (sum(rmaObject[[1]]$sigma2) + (rmaObject[[1]]$k - rmaObject[[1]]$p)/sum(diag(P))), 2)
  
  studyID <- rmaObject[[1]]$mf.r[[1]]$paperID
  resultID <- rmaObject[[1]]$mf.r[[1]]$effectID
  resR <- rma.mv(yi = rmaObject[[1]]$yi, V = rmaObject[[1]]$vi, random = ~ 1|studyID/resultID)
  resF <- rma.mv(yi = rmaObject[[1]]$yi, V = rmaObject[[1]]$vi)
  
  # Jackson's approach to I^2
  JI2 <- round(c(100 * (vcov(resR)[1,1] - vcov(resF)[1,1])/vcov(resR)[1,1]), 2)
  
  # Intra-class correlation of underlying true effects
  icc <- round(rmaObject[[1]]$sigma2[1]/sum(rmaObject[[1]]$sigma2), 2)
  
  c("Tau" = tau,
    "I^2" = I2,
    "Jackson's I^2" = JI2,
    "Between-cluster heterogeneity" = BW.hetero[1],
    "Within-cluster heterogeneity" = BW.hetero[2],
    "ICC" = icc)
}

# Proportion of significant results ---------------------------------------

propSig <- function(p.values = NA){
  as.integer(table(p.values < .05)[2])/length(p.values < .05)
}

# Permutation p-curve -----------------------------------------------------

# Subseting only the effects that are focal for the published paperID
pcurvePerm <- function(data, esEstimate = FALSE, plot = FALSE, nIterations = nIterationsPcurve){
  resultIDpcurve <- list(NA)
  resultPcurve <- matrix(ncol = 11, nrow = nIterationsPcurve)
  set.seed(1)
  for(i in 1:nIterationsPcurve){
    datPcurve <- data
    metaPcurve <- tryCatch(metagen(TE = yi, seTE = sqrt(vi), data = datPcurve),
                           error = function(e) NULL)
    modelPcurve <- tryCatch(pcurveMod(metaPcurve, effect.estimation = esEstimate, plot = plot), 
                            error = function(e) NULL)
    if(is.null(modelPcurve)){
      next
    } else {
      resultPcurve[i,] <- c("iterationNo" = i, "rightskew" = modelPcurve$pcurveResults[1,], "flatness" = modelPcurve$pcurveResults[2,])  
    }
    resultIDpcurve[[i]] <- datPcurve$effectID
  }
    colnames(resultPcurve) <- c("iterationNo", "rightskew.pBinomial", "rightskew.zFull", "rightskew.pFull", "rightskew.zHalf", "rightskew.pHalf", "flatness.pBinomial", "flatness.zFull", "flatness.pFull", "flatness.zHalf", "flatness.pHalf")
    medianResultPcurve <- resultPcurve %>% data.frame() %>% na.omit() %>% arrange(rightskew.zFull) %>% slice(ceiling(n()/2)) %>% unlist()
    metaResultPcurve <- tryCatch(metagen(TE = yi, seTE = sqrt(vi), data = data[data$effectID %in% unlist(resultIDpcurve[medianResultPcurve["iterationNo"]]),]), error = function(e) NULL)
    metaResultPcurve <<- metaResultPcurve
    if(!is.null(metaResultPcurve)){
      pcurveMod(metaResultPcurve, effect.estimation = esEstimate, plot = plot) 
    }
}

# Multiple-parameter selection models -------------------------------------
# 4/3-parameter selection model (4PSM/3PSM)
selectionModel <- function(data, minNoPvals = minPvalues, nIteration = nIterations, fallback = FALSE, steps = c(.025, 1), deltas = NA){
  resultSM <- matrix(ncol = 8, nrow = nIteration)
  set.seed(1)
  for(i in 1:nIteration){
    mydat <<- data #mydat <<- data[!duplicated.random(data$paperID) & data$focal == 1,]
    res <- rma(yi, vi,  data = mydat, method = "ML")
    # if <= min.pvalues p-values in an interval: return NULL
    pTable <- table(cut(mydat$p, breaks = c(0, .05, 0.5, 1)))
    if(fallback == TRUE | any(pTable < minNoPvals) | !anyNA(deltas)){
      threeFit <- tryCatch(selmodel(res, type = "stepfun", steps = steps, delta = deltas, alternative = "greater"),
                           error = function(e) NULL)
      threeOut <- if(is.null(threeFit)){
        next
      } else {
        round(c("est" = threeFit$beta, "se" = threeFit$se, "zvalue" = threeFit$zval, "pvalue" = threeFit$pval, "ciLB" = threeFit$ci.lb, "ciUB" = threeFit$ci.ub, "k" = threeFit$k, "steps" = length(threeFit$steps)), 3)
      }
      out <- threeOut
    } else { 
      fourFit <- tryCatch(selmodel(res, type = "stepfun", steps = c(.025, .5, 1), alternative = "greater"),
                          error = function(e) NULL)
      fourOut <- c("est" = fourFit$beta, "se" = fourFit$se, "zvalue" = fourFit$zval, "pvalue" = fourFit$pval, "ciLB" = fourFit$ci.lb, "ciUB" = fourFit$ci.ub, "k" = fourFit$k, "steps" = length(fourFit$steps))  
      if (is.null(fourFit)){
        out <- threeOut
      } else {
        out <- fourOut %>% round(., 3)
      }
    }
    resultSM[i,] <- out  
  }
  colnames(resultSM) <- c("est", "se", "zvalue", "pvalue", "ciLB", "ciUB", "k", "steps")
  resultSM <- resultSM %>% data.frame() %>% na.omit() %>% arrange(est) %>% slice(ceiling(n()/2)) %>% unlist()
  resultSM <<- resultSM
  resultSM
}

# Vevea & Woods step function model using a priori defined selection weights
veveaWoodsSM <- function(data, stepsDelta, nIteration = nIterationVWsensitivity){
  set.seed(1)
  do.call(rbind, lapply(stepsDelta[-1], function(delta) suppressWarnings(selectionModel(data, steps = stepsDelta$steps, deltas = delta, nIteration = nIterationVWsensitivity))))
}

# Robust Bayesian meta-analysis
bma <- function(data, seedNo = 1, chainsNo = robmaChains, nIterationBMA = robmaSamples){
  tryCatch(summary(data %>% filter(!is.na(yi) & !is.na(vi)) %$% RoBMA(d = yi, v = vi, study_names = label, seed = seedNo,
                                                                      chains = chainsNo, sample = nIterationBMA, burnin = ifelse(2*nIterationBMA/5 < 50, 50, 2*nIterationBMA/5), parallel = TRUE)), 
           error = function(e) NULL)
}

# PET-PEESE ---------------------------------------------------------------

#PET-PEESE with 4/3PSM as the conditional estimator instead of PET. 
# Also implemented the modified sample-size based estimator (see https://www.jepusto.com/pet-peese-performance/).
petPeese <- function(data, nBased = FALSE, selModAsCondEst = condEst){  # if nBased = TRUE, use the sample-size-based estimator, if FALSE, use the ordinary SE/var. If selModAsCondEst = TRUE, use the selection model as conditional estimator, otherwise use PET.
  viMatrix <- data %$% impute_covariance_matrix(vi, cluster = paperID, r = rho)  # compute the covariance matrix for the CHE working model
  
  if(nBased == TRUE){
    pet <<- robust.rma.mv(rma.mv(yi = yi ~ sqrt(nTerm), V = viMatrix, random = ~ 1|paperID/effectID, method = rmaMethod, sparse = TRUE, data = data), cluster = data$paperID)
  } else {
    pet <<- robust.rma.mv(rma.mv(yi = yi ~ sqrt(vi), V = viMatrix, random = ~ 1|paperID/effectID, method = rmaMethod, sparse = TRUE, data = data), cluster = data$paperID)
  }
  pet.out <- round(c(pet$b[1], pet$se[1], pet$zval[1], pet$pval[1], pet$ci.lb[1], pet$ci.ub[1]), 3)
  names(pet.out) <- c("PET estimate", "se", "zvalue", "pvalue", "ciLB", "ciUB")
  
  if(nBased == TRUE){
    peese <<- robust.rma.mv(rma.mv(yi = yi ~ nTerm, V = viMatrix, random = ~ 1|paperID/effectID, method = rmaMethod, sparse = TRUE, data = data), cluster = data$paperID)
  } else {
    peese <<- robust.rma.mv(rma.mv(yi = yi ~ vi, V = viMatrix, random = ~ 1|paperID/effectID, method = rmaMethod, sparse = TRUE, data = data), cluster = data$paperID)
  } 
  peese.out <- round(c(peese$b[1], peese$se[1], peese$zval[1], peese$pval[1], peese$ci.lb[1], peese$ci.ub[1]), 3)
  names(peese.out) <- c("PEESE estimate", "se", "zvalue", "pvalue", "ciLB", "ciUB")
  
  if(selModAsCondEst == TRUE){
    ifelse(resultSM["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * resultSM["est"] > 0,
           return(c(peese.out, pet.out)),  return(c(pet.out, peese.out)))
  } else {
    ifelse(pet$pval[1] < alpha & ifelse(exists("side") & side == "left", -1, 1) * pet$b[1] > 0,
           return(c(peese.out, pet.out)),  return(c(pet.out, peese.out)))
  } 
}

# WAAP-WLS estimator ------------------------------------------------------

# Code by Felix Schonbrodt and Evan Carter; https://github.com/nicebread/meta-showdown/blob/master/MA-methods/6-WAAP.R
# based on Stanley, T. D., Doucouliagos, H., & Ioannidis, J. P. A. (2017). Finding the power to reduce publication bias. Statistics in Medicine, 54(3), 30–19. http://doi.org/10.1002/sim.7228

# WLS estimator

# Weighted least squares estimator: Stanley, T. D., & Doucouliagos, H. (2015). Neither fixed nor random: weighted least squares meta-analysis. Statistics in Medicine, 34(13), 2116–2127. http://doi.org/10.1002/sim.6481
WLS.est <- function(yi, vi, long = FALSE) {
  se <- sqrt(vi)
  yi.precision <- 1/se
  yi.stand <- yi/se
  l1 <- lm(yi.stand ~ 0 + yi.precision)
  s1 <- summary(l1)
  res <- data.frame(
    method = "WLS",
    term = "b0",
    estimate = 	s1$coefficients[1, 1],
    std.error = s1$coefficients[1, 2],
    statistic = s1$coefficients[1, 2],
    p.value = s1$coefficients[1, 4],
    conf.low = confint(l1)[1],
    conf.high = confint(l1)[2]
  )
  returnRes(res, long = FALSE)
}

# return object: type = 1: WAAP, type = 2: WLS (must be numeric, otherwise it distorts the structure of the results object)
waapWLS <- function(yi, vi, est = c("WAAP-WLS"), long = FALSE) {
  
  # 1. determine which studies are in the top-N set
  
  # FE model (or, equivalently, WLS) as the proxy for true effect.
  WLS.all  <- WLS.est(yi, vi, long=FALSE)
  true.effect <- WLS.all$estimate
  
  # only select studies that are adequatly powered (Stanley uses a two-sided test)
  powered <- true.effect/2.8 >= sqrt(vi)
  
  # 2. compute the unrestricted weighted average (WLS) rma of either all or only adequatly powered studies
  # Combined estimator: WAAP-WLS	
  kAdequate <- sum(powered,na.rm=T)
  
  if (kAdequate >= 2) {
    res <- WLS.est(yi[powered], vi[powered], long=FALSE)
    res$method <- "WAAP-WLS"
    res <- plyr::rbind.fill(res, data.frame(method="WAAP-WLS", term="estimator", type=1, kAdequate=kAdequate))
  } else {
    res <- WLS.all
    res$method <- "WAAP-WLS"
    res <- plyr::rbind.fill(res, data.frame(method="WAAP-WLS", term="estimator", type=2, kAdequate=kAdequate))
  }
  returnRes(res, long = FALSE)
}


# Publication bias summary function-------------------------------

bias <- function(data = NA, rmaObject = NA, runRobMA = FALSE){
  # Correlation between the ES and precision (SE)
  esPrec <- cor(rmaObject[[1]]$yi, sqrt(rmaObject[[1]]$vi), method = "kendall")
  
  # Small-paperID effects correction
  # 3-parameter selection model
  resultSelModel <- selectionModel(data, minNoPvals = minPvalues, nIteration = nIterations, fallback = fallback)
  
  # Vevea & Woods selection model
  resultsVeveaWoodsSM <- veveaWoodsSM(data, stepsDelta)
  
  # Robust Bayesian model-averaging approach
  bmaMod <- if(runRobMA == TRUE){bma(data)}
  
  # PET-PEESE
  petPeeseOut <- petPeese(data)
  
  # WAAP-WLS
  waapWLSout <- data %$% waapWLS(yi, vi)
  waapWLSout[1, 9:10] <- waapWLSout[2, 9:10]
  waapWLSout <- waapWLSout[1,]
  
  # Permutation p-curve
  pcurvePermOut <- pcurvePerm(data, esEstimate = FALSE, plot = FALSE)
  
  # p-uniform* (van Aert & van Assen, 2021)
  resultPuniform <- matrix(ncol = 4, nrow = nIterations)
  set.seed(1)
  for(i in 1:nIterations){
    modelPuniform <- data %$% puni_star(yi = yi, vi = vi, alpha = alpha, side = side, method = "ML")
    resultPuniform[i,] <- c("est" = modelPuniform[["est"]], "ciLB" = modelPuniform[["ci.lb"]], "ciUB" = modelPuniform[["ci.ub"]], "p-value" = modelPuniform[["pval.0"]])
  }
    colnames(resultPuniform) <- c("est", "ciLB", "ciUB", "pvalue")
    puniformOut <- resultPuniform %>% data.frame() %>% na.omit() %>% arrange(est) %>% slice(ceiling(n()/2)) %>% unlist()
    puniformOut <<- resultPuniform
    
  # Publication bias results
  return(list("ES-precision correlation" = esPrec,
              "4/3PSM" = resultSelModel,
              "Vevea & Woods SM" = resultsVeveaWoodsSM,
              "PET-PEESE" = petPeeseOut,
              "WAAP-WLS" = waapWLSout,
              "p-uniform*" = puniformOut,
              "p-curve" = pcurvePermOut))
}

# Summary results ---------------------------------------------------------

maResults <- function(rmaObject = NA, data = NA, bias = T){
  list(
    "RMA results with model-based SEs" = rmaObject[[2]],
    "RVE SEs with Satterthwaite small-sample correction" = list("test" = coef_test(rmaObject[[2]], vcov = "CR2", cluster = data$paperID), "CIs" = conf_int(rmaObject[[2]], vcov = "CR2", cluster = data$paperID)),
    "Prediction interval" = pi95(rmaObject),
    "Heterogeneity" = heterogeneity(rmaObject),
    "Proportion of significant results" = propSig(data$p),
    "Publication bias" = if(bias ==T) {bias(data, rmaObject)} else {paste("Publication bias corrections not carried out")})
}

biasResults <- function(rmaObject = NA, data = NA){
  list(
    "Publication bias" = bias(data, rmaObject))
}

maResultsTable <- function(maResultsObject, metaAnalysis = TRUE, bias = TRUE){
  if(bias == TRUE & metaAnalysis == TRUE){
    noquote(c(
      "k" = as.numeric(maResultsObject[[1]]$k.all),
      "g [95% CI]" = paste(round(as.numeric(maResultsObject[[2]]$test$beta), 2), " [", round(maResultsObject[[2]]$CIs$CI_L, 2), ", ", round(maResultsObject[[2]]$CIs$CI_U, 2), "]", sep = ""),
      "95% PI [LB, UB]" = paste("[", round(maResultsObject$`Prediction interval`[1], 2), ", ", round(maResultsObject$`Prediction interval`[2], 2), "]", sep = ""),
      "SE" = round(maResultsObject[[2]]$test$SE, 2),
      round(maResultsObject[[4]]["Tau"], 2),
      "I^2" = paste(round(maResultsObject[[4]]["I^2"], 0), "%", sep = ""),
      "% significant" = paste(round(maResultsObject[["Proportion of significant results"]]*100, 0), "%", sep = ""),
      "3PSM est [95% CI]" = paste(round(maResultsObject[[6]][["4/3PSM"]]["est"], 2), " [", round(maResultsObject[[6]][["4/3PSM"]]["ciLB"], 2), ", ", round(maResultsObject[[6]][["4/3PSM"]]["ciUB"], 2), "]", sep = ""),
      "3PSM" = round(maResultsObject[[6]][["4/3PSM"]]["pvalue"], 3),
      "V&W [moderate/severe/extreme]" = paste(round(as.numeric(maResultsObject[[6]][["Vevea & Woods SM"]][1,1]), 2),"/",round(as.numeric(maResultsObject[[6]][["Vevea & Woods SM"]][2,1]), 2), "/", round(as.numeric(maResultsObject[[6]][["Vevea & Woods SM"]][3,1]), 2), sep = ""),
      "PET-PEESE est [95% CI]" = paste(round(maResultsObject[[6]][["PET-PEESE"]][1], 2), " [", round(maResultsObject[[6]][["PET-PEESE"]][5], 2), ", ", round(maResultsObject[[6]][["PET-PEESE"]][6], 2), "]", sep = ""),
      "PET-PEESE" = round(maResultsObject[[6]][["PET-PEESE"]][4], 3)
    ))
  } else if (metaAnalysis == FALSE & bias == TRUE) {
    noquote(c(
      "% significant" = paste(round(maResultsObject[["Proportion of significant results"]]*100, 0), "%", sep = ""),
      "3PSM est [95% CI]" = paste(round(maResultsObject[[1]][["4/3PSM"]]["est"], 2), " [", round(maResultsObject[[1]][["4/3PSM"]]["ciLB"], 2), ", ", round(maResultsObject[[1]][["4/3PSM"]]["ciUB"], 2), "]", sep = ""),
      "3PSM" = round(maResultsObject[[1]][["4/3PSM"]]["pvalue"], 3),
      "V&W [moderate/severe/extreme]" = paste(round(as.numeric(maResultsObject[[6]][["Vevea & Woods SM"]][1,1]), 2),"/",round(as.numeric(maResultsObject[[6]][["Vevea & Woods SM"]][2,1]), 2), "/", round(as.numeric(maResultsObject[[6]][["Vevea & Woods SM"]][3,1]), 2), sep = ""),
      "PET-PEESE est [95% CI]" = paste(round(maResultsObject[[1]][["PET-PEESE"]][1], 2), " [", round(maResultsObject[[1]][["PET-PEESE"]][5], 2), ", ", round(maResultsObject[[1]][["PET-PEESE"]][6], 2), "]", sep = ""),
      "PET-PEESE" = round(maResultsObject[[1]][["PET-PEESE"]][4], 3)
    ))
  } else {
    noquote(c(
      "k" = as.numeric(maResultsObject[[1]]$k.all),
      "g [95% CI]" = paste(round(as.numeric(maResultsObject[[2]]$test$beta), 2), " [", round(maResultsObject[[2]]$CIs$CI_L, 2), ", ", round(maResultsObject[[2]]$CIs$CI_U, 2), "]", sep = ""),
      "95% PI [LB, UB]" = paste("[", maResultsObject$`Prediction interval`[1], ", ", maResultsObject$`Prediction interval`[2], "]", sep = ""),
      "SE" = round(maResultsObject[[2]]$test$SE, 2),
      round(maResultsObject[[4]]["Tau"], 2),
      "I^2" = paste(round(maResultsObject[[4]]["I^2"], 0), "%", sep = "")
    ))
  }
}

# Return format function
# Code adapted from Carter, E. C., Schönbrodt, F. D., Hilgard, J., & Gervais, W. (2018). Correcting for bias in psychology: A comparison of meta-analytic methods. Retrieved from https://osf.io/rf3ys/.
# https://github.com/nicebread/meta-showdown/blob/master/MA-methods/7-Selection%20Models.R

# Return a effectID data frame either in wide or long format
returnRes <- function(res, long = TRUE, reduce = TRUE) {
  if (is.null(res)) return(NULL)
  
  # convert all factor columns to characters
  res %>% mutate_if(is.factor, as.character) -> res
  
  if (long == FALSE) {
    # return wide format
    return(res)
  } else {
    # transform to long format
    longRes <- melt(res, id.vars=c("method", "term"))
    if (reduce==TRUE & nrow(res) > 1) {longRes <- longRes %>% filter(!is.na(value)) %>% arrange(method, term, variable)}
    return(longRes)
  }
}

# t-test for summary statistics -------------------------------------------

tTestSummary <- function(mean1, mean2, sd1, sd2, n1, n2, withinSS = FALSE)
{
  if(withinSS == FALSE){
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt((1/n1 + 1/n2) * ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2)/(n1 + n2 - 2)) 
    df <- n1 + n2 - 2
    t <- (mean1 - mean2)/se 
    out <- c(mean1 - mean2, se, t, 2*pt(-abs(t),df))    
    names(out) <- c("Difference in means", "SE", "t-statistic", "p-value")
    return(out)  
  } else if(withinSS == TRUE){
    se <- sqrt((sd1^2 + sd2^2 - 2*sd1*sd2*rmCor)/(n1 - 1)) 
    df <- n1 - 1
    t <- (mean1 - mean2)/se 
    out <- c(mean1 - mean2, se, t, 2*pt(-abs(t), df))    
    names(out) <- c("Difference in means", "SE", "t-statistic", "p-value")
    return(out)
  }
}

# Random selection of effects ---------------------------------------------

# Choose effects from a single paperID by random (for the purpose of permutation-based methods)
duplicated.random = function(x, incomparables = FALSE, ...)
{
  if (is.vector(x))
  {
    permutation = sample(length(x))
    x.perm      = x[permutation]
    effectID.perm = duplicated(x.perm, incomparables, ...)
    effectID      = effectID.perm[order(permutation)]
    return(effectID)
  }
  else (is.matrix(x))
  {
    permutation = sample(nrow(x))
    x.perm      = x[permutation,]
    effectID.perm = duplicated(x.perm, incomparables, ...)
    effectID      = effectID.perm[order(permutation)]
    return(effectID)
  }
}

# GRIM & GRIMMER Output -----------------------------------------------------
grim <- function(dat){
  dat <- dat %>% mutate(items = ifelse(is.na(items), 0, items))
  outGrimM1 <- NA
  outGrimM2 <- NA
  datGRIM <- dat %>% filter(complete.cases(nExp, nCtrl, mExp, mCtrl, items))
  for(i in 1:nrow(datGRIM)){
    outGrimM1[i] <- grimTest(n = datGRIM[i,]$nExp, mean = datGRIM[i,]$mExp, items = datGRIM[i,]$items, decimals = 2)
    outGrimM2[i] <- grimTest(n = datGRIM[i,]$nCtrl, mean = datGRIM[i,]$mCtrl, items = datGRIM[i,]$items, decimals = 2)
  }
  
  datGRIM$outGrimM1 <- outGrimM1
  datGRIM$outGrimM2 <- outGrimM2
  datGRIM$inconsistenciesCountGRIM <- datGRIM %$% abs(outGrimM1 + outGrimM2 - 2)
  
  dat <<- datGRIM %>% 
    select(effectID, inconsistenciesCountGRIM) %>%
    left_join(dat, ., by = "effectID", keep = FALSE)
}

grimmer <- function(dat){
  dat <- dat %>% mutate(items = ifelse(is.na(items), 0, items))
  outGrimmerSD1 <- NA
  outGrimmerSD2 <- NA
  datGRIM <- dat %>% filter(complete.cases(nExp, nCtrl, mExp, mCtrl, sdExp, sdCtrl, items))
  for(i in 1:nrow(datGRIM)){
    outGrimmerSD1[i] <- grimmerTest(n = datGRIM[i,]$nExp, mean = datGRIM[i,]$mExp, SD = datGRIM[i,]$sdExp, items = datGRIM[i,]$items, decimals_mean = 2, decimals_SD = 2)
    outGrimmerSD2[i] <- grimmerTest(n = datGRIM[i,]$nCtrl, mean = datGRIM[i,]$mCtrl, SD = datGRIM[i,]$sdCtrl, items = datGRIM[i,]$items, decimals_mean = 2, decimals_SD = 2)
  }
  datGRIM$outGrimmerSD1 <- outGrimmerSD1
  datGRIM$outGrimmerSD2 <- outGrimmerSD2
  datGRIM$inconsistenciesCountGRIMMER <- datGRIM %$% ifelse((is.na(seExp) & is.na(seCtrl)), abs(outGrimmerSD1 + outGrimmerSD2 - 2), NA)
  
  dat <<- datGRIM %>% 
    select(effectID, inconsistenciesCountGRIMMER) %>%
    left_join(dat, ., by = "effectID", keep = FALSE)
}
# General Grim Test -------------------------------------------------------

# Code adapted from https://osf.io/scpbz/ , by Nick Brown and 
# https://aurelienallard.netlify.com/post/anaytic-grimmer-possibility-standard-deviations/, by Aurélien Allard

grimTest <- function (n, mean, items = 1, decimals = 2) {
  # if(n>10^decimals){
  #   print("The sample size is too big compared to the precision of the reported mean, it is not possible to apply GRIM.")
  # } else {
  if(items == 0 | is.na(items)){
    return(NA)} else {
      N <- n*items
      dust <- 1e-12
      gMean <- mean
      int <- round(mean * N) # nearest integer; doesn't matter if this rounds up or down
      frac <- int / N
      dif <- abs(mean - frac)
      gran <- ((0.1 ^ decimals) / 2) + dust # allow for rounding errors
      gMean <- round(int / N, decimals)
      consistent <- ifelse(gMean == mean, TRUE, FALSE)
      return(consistent)
      #}
    }
}

# General Grimmer Test ----------------------------------------------------

# Original GRIMMER effectID: -1 = GRIM inconsistent, 0 = GRIMMER inconsistent, 1 = mean & sd consistent
# For the present purposes, changed so that GRIM inconsistent returns 0 (GRIMMER inconsistent)
grimmerTest <- function(n, mean, SD, items = 1, decimals_mean = 2, decimals_SD = 2){
  # 
  # if(n>10^decimals_mean){
  #   print("Reported precision of mean too low, given N")
  # } else {
  # 
  #Applies the GRIM test, and computes the possible mean.
  if(items == 0 | is.na(items)){
    return(NA)} else {
      N <- n*items
      sum <- mean*N
      realsum <- round(sum)
      realmean <- realsum/N
      
      # Creates functions to round a number consistently up or down, when the last digit is 5
      
      round_down <- function(number, decimals=2){
        is_five <- number*10^(decimals+1)-floor(number*10^(decimals))*10
        number_rounded <- ifelse(is_five==5, floor(number*10^decimals)/10^decimals, round(number, digits = decimals))
        return(number_rounded)
      }
      
      round_up <- function(number, decimals=2){
        is_five <- number*10^(decimals+1)-floor(number*10^(decimals))*10
        number_rounded <- ifelse(is_five==5, ceiling(number*10^decimals)/10^decimals, round(number, digits = decimals))
        return(number_rounded)
      }
      
      # Applies the GRIM test, to see whether the reconstituted mean is the same as the reported mean (with both down and up rounding)
      
      consistency_down <- round_down(number = realmean, decimals = decimals_mean)==mean
      consistency_up <- round_up(number = realmean, decimals = decimals_mean)==mean
      
      if(consistency_down+consistency_up==0){
        return(0)
      }
      
      #Computes the lower and upper bounds for the sd.
      
      Lsigma <- ifelse(SD<5/(10^decimals_SD), 0, SD-5/(10^decimals_SD))
      Usigma <- SD+5/(10^decimals_SD)
      
      #Computes the lower and upper bounds for the sum of squares of items.
      
      Lowerbound <- (N-1)*Lsigma^2+N*realmean^2
      Upperbound <- (N-1)*Usigma^2+N*realmean^2
      
      #Checks that there is at least an integer between the lower and upperbound
      
      FirstTest <- ifelse(ceiling(Lowerbound)>floor(Upperbound), FALSE, TRUE)
      
      if(FirstTest==FALSE){
        return(0)
      }
      
      #Takes a vector of all the integers between the lowerbound and upperbound
      
      Possible_Integers <- ceiling(Lowerbound):floor(Upperbound)
      
      #Creates the predicted variance and sd
      
      Predicted_Variance <- (Possible_Integers-N*realmean^2)/(N-1)
      Predicted_SD <- sqrt(Predicted_Variance)
      
      #Computes whether one Predicted_SD matches the SD (trying to round both down and up)
      
      Rounded_SD_down <- round_down(Predicted_SD, decimals_SD)
      Rounded_SD_up <- round_up(Predicted_SD, decimals_SD)
      
      Matches_SD <- Rounded_SD_down==SD | Rounded_SD_up==SD
      
      if(sum(Matches_SD)==0){
        return(0)
      }
      
      #Computes first whether there is any integer between lower and upper bound, and then whether there is 
      #an integer of the correct oddness between the lower and upper bounds.
      oddness <- realsum%%2
      Matches_Oddness <- Possible_Integers%%2==oddness
      Third_Test <- Matches_SD&Matches_Oddness
      return(ifelse(
        sum(Third_Test)==0, 0, 1)
      )
    }
}

