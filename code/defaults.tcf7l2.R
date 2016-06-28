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
    transformation = NULL,
    covariates = c("sex","bw3","group4B"),
    outlier.function = function (x) x > 60,
    outliers = c(49978, 50933, 50986, 51063, 51093, 52403, 52428, 53135,
                 53161, 53268, 53471, 53476, 53513, 53517, 53577, 53631),
    gemma = c("all", "wt", "het")),
  
  fastglucose = list(
    transformation = NULL,
    covariates = c("sex","bw3"),
    outlier.function = function (x) x > 60,
    outliers = c(50964, 51031, 51032, 53113, 53146, 53177,
                 53241, 53532, 53633, 53647),
    gemma = c("all", "wt", "het")),
  
  bw2 = list(
    transformation = NULL,
    covariates = c("sex","littersize"),
    outlier.function = NULL,
    outliers = c(50901,53186,53517),
    gemma = c("all", "wt", "het")),
  
  bw3 = list(
    transformation = NULL,
    covariates = c("sex", "littersize"),
    outlier.function = NULL,
    outliers = c(50906, 51090, 52417, 53102, 53146, 53148, 53151, 53256, 53422,
                 53487, 53517, 53595),
    gemma = c("all", "wt", "het")),
  
  weightgain = list(
    transformation = NULL,
    covariates = c("sex"),
    outlier.function = function (x) x > 13,
    outliers = c(50906, 51059, 51060, 53102, 53132, 53151, 53422, 53517,
                 53595, 53612, 53613),
    gemma = c("all", "wt", "het")),
  
  d2context = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("sex", "FCtimeofday", "agouti", "d1tone"),
    outlier.function = function (x) x > 2,
    outliers = c(50970, 53104, 53141, 53652),
    gemma = c("all", "wt", "het")),
  
  d3altcontext = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("sex", "FCtimeofday", "agouti"),
    outlier.function = NULL,
    outliers = c(50964, 50970, 50973, 51078, 51082, 52439,
                 53472, 53509, 53568, 53629, 53697),
    gemma = c("all", "wt", "het")),
  
  d3tone = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("sex", "FCtimeofday", "agouti"),
    outlier.function = NULL,
    outliers = c(50973, 51021, 52409, 52415, 53110, 53426, 53486,
                 53509, 53558, 53613, 53697),
    gemma = c("all", "wt", "het")),
   
  PPI3 = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("PPIbox3"),
    outlier.function = function (x) x < (-1),
    outliers = c(49980, 49980, 51015, 51015, 51071, 51075, 51090, 51090, 
                 52437, 52437, 53102, 53174, 53260, 53263, 53263, 53552),
    gemma = c("all", "wt", "het")),
  
  PPI6 = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("PPIbox3"),
    outlier.function = function (x) x < (-1),
    outliers = c(49980, 51015, 51078, 52444, 53487, 53537, 53594,
                 53647, 53649),
    gemma = c("all", "wt", "het")),
  
  PPI12 = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("PPIbox3"),
    outlier.function = function (x) x < (-0.95),
    outliers = c(49980, 51056, 53291, 53649),
    gemma = c("all", "wt", "het")),
  
  PPIavg = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("PPIbox3"),
    outlier.function = function (x) x < (-1),
    outliers = c(50967, 51015, 51034, 52417, 52437, 53123, 53126, 53504,
                 53541, 53647, 53649, 53667),
    gemma = c("all", "wt", "het")),
  
  PPIstartle = list(
    transformation = function (x) log10(x),
    covariates = c("sex", "bw2", "PPIbox1", "PPIbox2", "PPIbox3", "PPIbox4"),
    outlier.function = NULL,
    outliers = c(50927, 51015, 52406, 53181, 53250, 53278, 53431, 53487),
    gemma = c("all", "wt", "het")),
  
  PPIhabit = list(
    transformation = function(x) log10(x),
    covariates = c("p120b1", "sex", "bw2", "PPIbox1", "PPIbox2", "PPIbox3", "PPIbox4"),
    outlier.function = NULL,
    outliers = c(50938, 51040, 51086, 52427, 53129, 53208, 53618, 53674),
    gemma = c("all", "wt", "het")),
  
  centerduration = list(
    transformation = function (x) logit10(project.onto.interval(x / 400, 0.001, 0.999)),
    covariates = c("sex"),
    outlier.function = function (x) x < (-1.4) | x > 1.5,
    outliers = c(50985, 51019, 51021, 52446, 53196, 53217, 53420,
                 53429, 53512, 53570),
    gemma = c("all", "wt", "het")),
  
  verticalactivity = list(
    # Note: this transformation depends on totalactivity data from prepared.pheno
    transformation = function(x) {
        logit10(project.onto.interval(x / prepared.pheno$totalactivity, 0.001, 0.999))
    },
    covariates = c("sex", "bw2", "oftbox9", "oftbox10", "oftbox11"),
    outlier.function = function (x) x < (-0.75) | x > 1,
    outliers = c(50925, 50992, 51026, 51091, 53177, 53242, 53245, 53257,
                 53294, 53430, 53490, 53523),
    gemma = c("all", "wt", "het")),
  
  totalactivity = list(
    transformation = function(x) log10(x),
    covariates = NULL,
    outlier.function = function (x) x < (-0.5) | x > (0.5),
    outliers = c(51003, 53114, 53163, 53185, 53202, 53522, 53653),
    gemma = c("all", "wt", "het")),
  
  fstdist = list(
    transformation = function(x) log10(x),
    covariates = c("sex", "agouti"),
    outlier.function = function (x) x < (-0.2) | x > 0.2,
    outliers = c(50959, 50961, 51085, 53680),
    gemma = c("all", "wt", "het")),
  
  immobfreq = list(
    transformation = NULL,
    covariates = c("sex", "bw2"),
    outlier.function = function (x) x > 400,
    outliers = c(53437, 51034, 53662, 50994),
    gemma = c("all", "wt", "het")),
  
  absimmobdur = list(
    transformation = function(x) log10(pmax(0.01, x)),
    covariates = c("agouti"),
    outlier.function = function (x) x < (-2),
    outliers = c(50983, 51001, 53115, 53401, 53661),
    gemma = c("all", "wt", "het")),
  
  absimmoblatency = list(
    transformation = function(x) log10(pmax(0.01, x)),
    covariates = NULL,
    outlier.function = function (x) x < (-1),
    outliers = c(51006, 53115, 53226, 53294, 53296, 53504, 53586, 53661),
    gemma = c("all", "wt", "het"))

)



