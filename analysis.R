#' ---
#' title: "SUPERAR - reanalysis of meta-analyses on the effect of music education"
#' author: "Ivan Ropovik"
#' output:
#'    html_document:
#'       code_folding: hide
#'       toc: true
#'       toc_float: true
#'       fig_retina: 2
#'       theme: paper
#' always_allow_html: yes
#' ---
knitr::opts_chunk$set(echo = FALSE, warning = FALSE)

#' **This is the supplementary analytic output for the reanalysis of meta-analyses on the effect of music education, carried out for the SUPERAR project**
#' 
#' ------------------------------------
#' 
#' **Brief information about the methods used in the analysis:**
#' 
#' **RMA results with model-based SEs**
#' k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
#'
#' **RVE SEs with Satterthwaite small-sample correction**
#' Estimate based on a multilevel RE model with constant sampling correlation model (CHE - correlated hierarchical effects - working model) (Pustejovsky & Tipton, 2020; https://osf.io/preprints/metaarxiv/vyfcj/). 
#' Interpretation of naive-meta-analysis should be based on these estimates.
#'
#' **Prediction interval**
#' Shows the expected range of true effects in similar studies.
#' As an approximation, in 95% of cases the true effect in a new *published* study can be expected to fall between PI LB and PI UB.
#' Note that these are non-adjusted estimates. An unbiased newly conducted study will more likely fall in an interval centered around bias-adjusted estimate with a wider CI width.
#'
#' **Heterogeneity**
#' Tau can be interpreted as the total amount of heterogeneity in the true effects. 
#' I^2$ represents the ratio of true heterogeneity to total variance across the observed effect estimates. Estimates calculated by two approaches are reported.
#' This is followed by separate estimates of between- and within-cluster heterogeneity and estimated intra-class correlation of underlying true effects.
#' 
#' **Proportion of significant results**
#' What proportion of effects were statistically at the alpha level of .05.
#' 
#' **ES-precision correlation**
#' Kendalls's correlation between the ES and precision.
#' 
#' **4/3PSM**
#' Applies a permutation-based, step-function 4-parameter selection model (one-tailed p-value steps = c(.025, .5, 1)). 
#' Falls back to 3-parameter selection model if at least one of the three p-value intervals contains less than 5 p-values.
#' For this meta-analysis, we applied 3-parameter selection model by default as there were only 11 independent effects in the opposite direction overall (6%), causing the estimates to be unstable across iterations.
#' pvalue = p-value testing H0 that the effect is zero. ciLB and ciUB are lower and upper bound of the CI. k = number of studies. steps = 3 means that the 4PSM was applied, 2 means that the 3PSM was applied.
#' We also ran two sensitivity analyses of the selection model, the Vevea & Woods (2005) step function model with a priori defined selection weights and the Robust Bayesian Meta-analysis model employing the model-averaging approach (Bartoš & Maier, 2020).
#' 
#' **PET-PEESE**
#' Estimated effect size of an infinitely precise study. Using 4/3PSM as the conditional estimator instead of PET (can be changed to PET). If the PET-PEESE estimate is in the opposite direction, the effect can be regarded nil. 
#' By default (can be changed to PET), the function employs a modified sample-size based estimator (see https://www.jepusto.com/pet-peese-performance/). 
#' It also uses the same RVE sandwich-type based estimator in a CHE (correlated hierarchical effects) working model with the identical random effects structure as the primary (naive) meta-analytic model. 
#' 
#' We report results for both, PET and PEESE, with the first reported one being the primary (based on the conditional estimator).
#' 
#' **WAAP-WLS**
#' The combined WAAP-WLS estimator (weighted average of the adequately powered - weighted least squares) tries to identify studies that are adequately powered to detect the meta-analytic effect. 
#' If there is less than two such studies, the method falls back to the WLS estimator (Stanley & Doucouliagos, 2015). If there are at least two adequately powered studies, WAAP returns a WLS estimate based on effects from only those studies.
#' 
#' type = 1: WAAP estimate, 2: WLS estimate. kAdequate = number of adequately powered studies
#' 
#' **p-uniform**
#' P-uniform* is a selection model conceptually similar to p-curve. It makes use of the fact that p-values follow a uniform distribution at the true effect size while it includes also nonsignificant effect sizes.
#' Permutation-based version of p-uniform method, the so-called p-uniform* (van Aert, van Assen, 2021).
#' 
#' **p-curve**
#' Permutation-based p-curve method. Output should be self-explanatory. For more info see p-curve.com
#' 
#' **Power for detecting SESOI and bias-corrected parameter estimates**
#' Estimates of the statistical power for detecting a smallest effect sizes of interest equal to .20, .50, and .70 in SD units (Cohen's d). 
#' A sort of a thought experiment, we also assumed that population true values equal the bias-corrected estimates (4/3PSM or PET-PEESE) and computed power for those.
#' 
#' **Handling of dependencies in bias-correction methods**
#' To handle dependencies among the effects, the 4PSM, p-curve, p-uniform are implemented using a permutation-based procedure, randomly selecting only one focal effect (i.e., excluding those which were not coded as being focal) from a single study and iterating nIterations times.
#' Lastly, the procedure selects the result with the median value of the ES estimate (4PSM, p-uniform) or median z-score of the full p-curve (p-curve).


# Settings ----------------------------------------------------------------
#+ setup, include = FALSE
rm(list = ls())
rmaMethod <- "REML"

# Assumed constant sampling correlation
rho <- 0.5

# For a more elaborate output from the pub bias tests/correction methods, set briefBias to FALSE
biasOn <- F
briefBias <- TRUE

# Side argument for the p-uniform* and conditional estimator of PET-PEESE. If the target effect should be in negative values, set to "left", otherwise "right".
side <- "right"

# Define whether to use one-tailed or two-tailed test for PET-PEESE, 3PSM, and p-uniform*.
# Recommended by Stanley (2016) for literature where small sample-size studies are rather the norm.
# Assuming alpha level of .05 for the two-tailed test
test <- "one-tailed"

# No of simulations for the permutation-based bias correction models and p-curve specifically
nIterations <- 1000 # Set to 5 just to make code checking/running fast.
nIterationsPcurve <- 20
nIterationVWsensitivity <- 20 # Number of iterations for the Vevea & Woods (2005) step function model sensitivity analysis 

# Number of chains and iterations for Robust Bayesian model-averaging approach
runRobMA <- FALSE
robmaChains <- 2
robmaSamples <- 1000

# Controls for PET-PEESE
condEst <- FALSE

# Controls for the multiple-parameter selection models 

# Whether to apply a 4- or 3-parameter selection model. If fallback == TRUE, the procedure falls back to the 3-parameter selection model. 
# This should be selected when too few effects in the opposite side make the estimate unstable.
fallback <- TRUE
# Even when fallback == FALSE, the 4-parameter selection model still falls back to 3 parameters for the given iteration if,
# (1) it fails to converge or (2) the number of p-values in each of the step intervals gets smaller than minPvalues.
minPvalues <- 3

# Steps and delta parameters for Vevea & Woods selection models 
# Can be adjusted if a different selection process is assumed. 
# Please note that steps vector represents one-tailed p-values.
stepsDelta <- data.frame(
  steps =     c(.0025, .005, .0125, .025, .05, .10, .25, .50, 1),
  moderateSelection = c(1, 0.99, 0.97, 0.95, 0.80, 0.60, 0.50, 0.50, 0.50),
  severeSelection =   c(1, 0.99, 0.97, 0.95, 0.65, 0.40, 0.25, 0.25, 0.25),
  extremeSelection =  c(1, 0.98, 0.95, 0.90, 0.50, 0.20, 0.10, 0.10, 0.10))

# Libraries
library(tidyr)
library(metafor)
library(compute.es)
library(psych)

# Source funcions
source("pcurvePlotOption.R")
source("functions.R")

# Define the path to the 'data' folder within the current working directory
dataFolderPath <- file.path(getwd(), "data")
excelFiles <- list.files(path = dataFolderPath, pattern = "\\.xlsx$", full.names = TRUE)

# Function to read all sheets from a given Excel file and ignore sheets with less than 5 rows
readExcelSheets <- function(fileName) {
  sheets <- excel_sheets(fileName)
  sheetsData <- lapply(sheets, function(sheet) {
    dataset <- read_excel(fileName, sheet = sheet)
    if (nrow(dataset) >= 2) {
      return(dataset)
    } else {
      return(NULL)
    }
  })
  # Filter out NULL entries (sheets with less than 5 rows) and corresponding sheet names
  validSheets <- !sapply(sheetsData, is.null)
  sheetsData <- sheetsData[validSheets]
  sheetNames <- sheets[validSheets]
  
  setNames(sheetsData, sheetNames)
}

# Create the nested list, excluding datasets with less than 5 rows
dat <- setNames(lapply(excelFiles, readExcelSheets), gsub("\\.xlsx$", "", basename(excelFiles)))

# Function to calculate variance and p-value for each dataset
calculateStats <- function(dataset) {
  # Add yi to the dataset if not present
  if (!("yi" %in% names(dataset))) {
    dataset$yi <- NA
  }
  
  # Convert yi to numeric if necessary
  dataset$yi <- as.numeric(dataset$yi)
  
  if ("vi" %in% names(dataset) && any(!is.na(dataset$yi))) {
    # If vi is in the dataset and any yi value is not NA, calculate p
    dataset <- dataset %>%
      mutate(
        p = ifelse(!is.na(vi), 2 * pnorm(-abs(yi / sqrt(vi))), NA) # Two-tailed p-value
      )
  }
  
  if (all(c("ciLB", "ciUB") %in% names(dataset)) && any(!is.na(dataset$yi))) {
    # Convert to numeric if necessary
    dataset$ciLB <- as.numeric(dataset$ciLB)
    dataset$ciUB <- as.numeric(dataset$ciUB)
    
    # Check if conversion was successful
    if (all(!is.na(dataset$ciLB), !is.na(dataset$ciUB))) {
      # Calculate variance and p-value if ciLB, ciUB, and yi are present
      dataset <- dataset %>%
        mutate(
          se = (ciUB - ciLB) / (2 * 1.96), # Standard error from CI
          vi = se^2, # Variance
          p = 2 * pnorm(-abs(yi / se)) # Two-tailed p-value
        )
    }
  }
  return(dataset)
}

# Apply the function to each dataset in the 'dat'
dat <- lapply(dat, function(file) {
  lapply(file, calculateStats)
})

# Function to apply rmaCustom, maResults, and maResultsTable to a dataset
analyzeDatasets <- function(dataset) {
  rmaOverall <- dataset %>% rmaCustom()
  resultsOverall <- dataset %>% maResults(rmaObject = rmaOverall, bias = biasOn)
  resultsTable <- maResultsTable(resultsOverall, bias = biasOn)
  list(resultsOverall = resultsOverall, resultsTable = resultsTable)
}

# Function to extract only the rmaOverall object
extractRmaObjects <- function(dataset) {
  rmaOverall <- dataset %>% rmaCustom()
  return(rmaOverall)
}

# Apply the analyzeDatasets function to each dataset in the 'dat' nested list
results <- lapply(dat, function(file) {
  lapply(file, analyzeDatasets)
})

# Apply the extractRmaOverall function to each dataset in the 'dat' nested list
rmaObjects <- lapply(dat, function(file) {
  lapply(file, function(dataset) {
    extractRmaObjectsResult <- extractRmaObjects(dataset)
    extractRmaObjectsResult[[1]]  # Extract only the first object
  })
})

#+ include = TRUE
#'# Results
results

#'## Forest plots
for (fileName in names(rmaObjects)) {
  for (sheetName in names(rmaObjects[[fileName]])) {
    plotTitle <- paste("Forest Plot for", fileName, "-", sheetName)
    metafor::forest(rmaObjects[[fileName]][[sheetName]], main = plotTitle)
  }
}

#'## Funnel plots
for (fileName in names(rmaObjects)) {
  for (sheetName in names(rmaObjects[[fileName]])) {
    plotTitle <- paste("Funnel Plot for", fileName, "-", sheetName)
    metafor::funnel(rmaObjects[[fileName]][[sheetName]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", main = plotTitle)
  }
}