# This file contains a list which specifies the default analysis
# settings for each phenotype recorded in the Tcf7l2 cohort:
#
#   (a) tranformation function
#   (b) covariate names
#   (c) outlier removal function
#   (d) specified outliers
#
model.info.tcf7l2 <- list(
  
  baseglucose = list(
    transformation   = NULL,
    covariates       = c("sex","bw3","group4B"),
    outlier.function = function (x) x > 60,
    outliers         = c(49978,50933,50986,51063,51093,52403,52428,53135,
                         53161,53268,53471,53476,53513,53517,53577,53631)),
  
  fastglucose = list(
    transformation   = NULL,
    covariates       = c("sex","bw3"),
    outlier.function = function (x) x > 60,
    outliers         = c(50964,51031,51032,53113,53146,53177,
                         53241,53532,53633,53647)),
  
  bw2 = list(
    transformation   = NULL,
    covariates       = c("sex","littersize"),
    outlier.function = NULL,
    outliers         = c(50901,53186,53517)),
  
  d2context = list(
    transformation   = function (x) logit10(project.onto.interval(x,1,99)/100),
    covariates       = c("sex","FCtimeofday","agouti","d1tone"),
    outlier.function = function (x) x > 2,
    outliers         = c(50970,53104,53141,53652)),
  
  d3tone = list(
    transformation   = function (x) logit10(project.onto.interval(x,1,99)/100),
    covariates       = c("sex","FCtimeofday","agouti"),
    outlier.function = NULL,
    outliers         = c(50973,51021,52409,52415,53110,53426,53486,
                         53509,53558,53613,53697)),
   
  PPI6 = list(
    transformation   = function (x) logit10(project.onto.interval(x,1,99)/100),
    covariates       = "PPIbox3",
    outlier.function = function (x) x < (-1),
    outliers         = c(49980,51015,51078,52444,53487,53537,53594,
                         53647,53649)),
  
  PPIstartle = list(
    transformation   = function (x) log10(x),
    covariates       = c("sex","bw2","PPIbox1","PPIbox2","PPIbox3","PPIbox4"),
    outlier.function = NULL,
    outliers         = c(50927,51015,52406,53181,53250,53278,53431,53487)),
  
  centerduration = list(
    transformation   =
      function (x) logit10(project.onto.interval(x/400,0.001,0.999)),
    covariates       = "sex",
    outlier.function = function (x) x < (-1.4) | x > 1.5,
    outliers         = c(50985,51019,51021,52446,53196,53217,53420,
                         53429,53512,53570)),
  
  totalactivity = list(
    transformation   = function(x) log10(x),
    covariates       = NULL,
    outlier.function = function (x) x < (-0.5) | x > (0.5),
    outliers         = c(51003,53114,53163,53185,53202,53522,53653)))
