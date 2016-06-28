# This file contains a list which holds the default settings for each phenotype's:
#   (a) tranformation function
#   (b) covariate names
#   (c) outlier removal function
#   (d) specified outliers
#   (e) which types of gemma to run

model.info <- list(

  d50bw = list(
    transformation = NULL,
    covariates = c("sex", "littersize"),
    outlier.function = NULL,
    outliers = c("55386", "55331", "55561", "55563", "57039", "55247", "55650", "55679",
                 "57036", "55063", "55347", "55317", "55832", "55087"),
    gemma = c("all", "wt", "het")),
  
  d100bw = list(
    transformation = NULL,
    covariates = c("sex", "littersize"),
    outlier.function = NULL,
    outliers = c("55658", "55657", "55842", "55028", "55027", "55331", "55665",
                 "55555", "55590", "57001", "55897", "55679", "55063", "55317",
                 "55320", "55832", "55087", "55125", "55057"),
    gemma = c("all", "wt", "het")),

  weightgain = list(
    transformation = NULL,
    covariates = c("sex"),
    outlier.function = NULL,
    outliers = c("55657", "55842", "55398", "55324", "55027", "55759", "55261",
                 "55193", "57001", "55897", "57036", "55112", "55031", "55156", "55057"),
    gemma = c("all", "wt", "het")),
    
  ldtotalcm = list(
    transformation = NULL,
    covariates = NULL,
    outlier.function = NULL,
    outliers = c("55615", "55289", "55531", "55039", "55329", "55761", "55180",
                 "55693", "55804", "55827", "55604", "55624", "55307", "55525"),
    gemma = c("all", "wt", "het")),  
  
  pctdurlight = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = NULL,
    outlier.function = NULL,
    outliers = c("55615", "55531", "55039", "55488", "55761", "55592", "55252",
                 "55804", "55874", "55009", "55827", "55087", "55080", "55604",
                 "55624", "55541", "55056", "55058"),
    gemma = c("all", "wt", "het")),

  latencytohidden = list(
    transformation = function (x) log10(x + 0.001),
    covariates = NULL,
    outlier.function = NULL,
    outliers = c("57063", "57049", "55361", "55359", "55709", "55379", "55042",
                 "55252", "55180", "55874", "55209", "55862", "55815", "55818"),
    gemma = c("all", "wt", "het")),

  d1totalactivity = list(
    transformation = function (x) log10(x),
    covariates = NULL,
    outlier.function = NULL,
    outliers = c("55614", "55359", "55799", "55570", "55689", "57024", "55234", "55246",
                 "55248", "55884", "55873", "55551", "55525"),
    gemma = c("all", "wt", "het")),

  d1centerdur = list(
    transformation = function (x) logit10(project.onto.interval(x / 400, 0.001, 0.999)),
    covariates = c("sex"),
    outlier.function = NULL,
    outliers = c("55109", "55359", "55405", "55323", "55595", "57056", "55826", "55269",
                 "57029", "55006", "55008", "55676", "55827", "55021", "55638"),
    gemma = c("all", "wt", "het")),

  d2totalactivity = list(
    transformation = function (x) logit10(
                                  project.onto.interval(x / 15300, .001, .999)),
    covariates = NULL,
    outlier.function = NULL,
    outliers = c("55514", "55577", "55570", "57024", "55223", "55893", "55318",
                 "55829", "55827", "55638", "55172", "55525", "55522"),
    gemma = c("all", "wt", "het")),

  d2centerdur = list(
    transformation = function (x) logit10(
                                  project.onto.interval(x / 454, 0.001, .999)),
    covariates = c("sex"),
    outlier.function = NULL,
    outliers = c("55359", "56133", "55297", "55793", "55568", "55562", "55248", "55246",
                 "55223", "55826", "55893", "55861", "55020", "55021", "55095", "55155"),
    gemma = c("all", "wt", "het")),
  
  d3.d2totalactivity = list(
    transformation = NULL,
    covariates = "d2totalactivity",
    outlier.function = NULL,
    outliers = c("55208", "55209", "57037", "55103", "55113", "55577", "55602", "55704", "55781", "55805", "57044", "55569", "55856", "55862"),
    gemma = c("all", "wt", "het")),
  
  fstdist = list(
    transformation = function (x) log10(x),
    covariates = c("sex", "agouti"),
    outlier.function = NULL,
    outliers = c("55811", "55385", "55463", "55370", "55673", "55846", "55884", "55825",
                 "55684", "55679", "55347", "55437", "55098"),
    gemma = c("all", "wt", "het")),
    
  highmobdur = list(
    transformation = function (x) log10(x + .001),
    covariates = c("agouti", "d100bw"),
    outlier.function = NULL,
    outliers = c("55290", "55482", "55028", "55328", "55792", "55846", "55824", "57028", 
                 "55203", "55001", "55875", "55770", "55862", "55437", "55006"),
    gemma = c("all", "wt", "het")),
  
  immobdur = list(
    transformation = function (x) log10(pmax(0.01, x)),
    covariates = c("agouti"),
    outlier.function = NULL,
    outliers = c("57063", "55663", "57062", "57066", "57067", "57068", "55110", "55107",
                 "55342", "55615", "55359", "55038", "55385", "55723", "55364", "55370", 
                 "55857", "55569", "55232", "55180", "55220", "55825", "55804", "55184",
                 "55186", "55100", "55137", "55098", "55631", "55503", "55045", "55056"),
    gemma = c("all", "wt", "het")),
    
  immobfreq = list(
    transformation = function (x) logit10(project.onto.interval(x / 139, .001, .999)),
    covariates = NULL,
    outlier.function = NULL,
    outliers = c("55615", "55359", "55482", "55038", "55385", "55673", "55857", "55569", 
                 "55825", "55804", "55184", "55186", "55100", "55098", "55631", "55503",
                 "55045", "55056"),
    gemma = c("all", "wt", "het")),
  
  ppi3 = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("ppibox3"),
    outlier.function = NULL,
    outliers = c("55659", "55662", "55405", "55498", "55564", "55745", "57038", "55692", 
                 "55063", "55320", "55005", "55677", "55305", "55603", "55095", "55156", 
                 "55161", "55501", "55060"),
    gemma = c("all", "wt", "het")),
  
  ppi6 = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("ppibox3"),
    outlier.function = NULL,
    outliers = c("55405", "55498", "55043", "55325", "55363", "55564", "57038", "55182", "55692", 
                 "55693", "55750", "55114", "55136", "55080"),
    gemma = c("all", "wt", "het")),
  
  ppi12 = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("ppibox1", "ppibox3", "ppibox5"),
    outlier.function = NULL,
    outliers = c("55405", "55043", "55596", "55564", "55745", "57038", "55182", "55692",
                 "55114", "55003", "55033", "55348", "55080", "55060"),
    gemma = c("all", "wt", "het")),
  
  ppiavg = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("ppibox1", "ppibox3", "ppibox4", "ppibox5"),
    outlier.function = NULL,
    outliers = c("55662", "55405", "55498", "55043", "55325", "55370", "55369", "55669",
                 "55564", "55745", "57038", "55182", "55692", "55693", "55114", "55003",
                 "55080", "55136", "55060"),
    gemma = c("all", "wt", "het")),
  
  startle = list(
    transformation = function (x) log10(x),
    covariates = c("sex", "d50bw", "ppibox1", "ppibox2", "ppibox3"),
    outlier.function = NULL,
    outliers = c("55573", "55043", "55553", "55292", "55648", "55866", "55694", "55693", 
                 "55184", "57029"),
    gemma = c("all", "wt", "het")),
  
  p120b4 = list(
    transformation = function (x) log10(x),
    covariates = c("p120b1", "sex", "d50bw",
                   "ppibox1", "ppibox2", "ppibox3", "ppibox4"),
    outlier.function = NULL,
    outliers = c("55252", "55259", "55179", "57004", "55865", "55694", "55693", "55409",
                 "55437", "55816"),
    gemma = c("all", "wt", "het")),
    
#   Note: For this phenotype (avtoned1), it makes more sense to implement the
#   outlier removal function rather than specify all of the outliers by ID. There
#   is an obvious threshold to use (residual less than around -1.5). Mice with 
#   residuals this low correspond to an untransformed avtoned1 value of 0. The transformation
#   puts these values all at -2 (which is log10(0.01)). This then corresponds to residuals
#   all less than ~ -1.5.

  avtoned1 = list(
    transformation = function (x) log10(pmax(0.01, x)),
    covariates = c("pretraind1"),
    outlier.function = function (x) x < (-1.5),
    outliers = NULL,
    gemma = c("all", "wt", "het")),
  
  pretraind1 = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = NULL,
    outlier.function = NULL,
    outliers = c("55590"),
    gemma = c("all", "wt", "het")),

  avcontextd2 = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("sex", "fcbox", "avtoned1"),
    outlier.function = NULL,
    outliers = c("57067", "55763", "55355", "55666", "55755", "55726", "55648", "55252",
                 "55825", "55680", "55581", "55114", "55099", "55133", "55075", "55077",
                 "55020", "55462"),
    gemma = c("all", "wt", "het")),

  avaltcontextd3 = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("sex", "fcbox"),
    outlier.function = NULL,
    outliers = c("55766", "55027", "55297", "55845", "55792", "55590", "55558", "55690",
                 "57058", "57042", "55438", "55305", "55537", "55520", "55059"),
    gemma = c("all", "wt", "het")),

  avtoned3 = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("sex", "fctimeofday"),
    outlier.function = NULL,
    outliers = c("55766", "55027", "55297", "55845", "55792", "55590", "55558", "55690",
                 "57058", "57042", "55438", "55305", "55537", "55520", "55059"),
    gemma = c("all", "wt", "het"))

)
