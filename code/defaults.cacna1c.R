# This list specifies the default analysis settings for each phenotype
# analyzed in the Cacna1c cohort:
#
#   (a) transformation function
#   (b) covariate names
#   (c) specified outliers
#
model.info.cacna1c <- list(

  pctdurlight = list(
    transformation   = function(x) logit10(project.onto.interval(x,1,99)/100),
    covariates       = NULL,
    outliers         = c(55615,55531,55039,55488,55761,55592,55252,
                         55804,55874,55009,55827,55087,55080,55604,
                         55624,55541,55056,55058)),

  d1totalactivity = list(
    transformation   = function (x) log10(x),
    covariates       = NULL,
    outliers         = c(55614,55359,55799,55570,55689,57024,55234,55246,
                         55248,55884,55873,55551,55525)),

  d1centerdur = list(
    transformation   = function (x) logit10(project.onto.interval(x/400,0.001,
                                                                  0.999)),
    covariates       = "sex",
    outliers         = c(55109,55359,55405,55323,55595,57056,55826,55269,
                         57029,55006,55008,55676,55827,55021,55638)),

  d3.d2totalactivity = list(
    transformation   = NULL,
    covariates       = "d2totalactivity",
    outliers         = c(55208,55209,57037,55103,55113,55577,55602,55704,
                         55781,55805,57044,55569,55856,55862)),
  
  immobdur = list(
    transformation = function (x) log10(pmax(0.01,x)),
    covariates       = "agouti",
    outliers         = c(57063,55663,57062,57066,57067,57068,55110,55107,
                         55342,55615,55359,55038,55385,55723,55364,55370, 
                         55857,55569,55232,55180,55220,55825,55804,55184,
                         55186,55100,55137,55098,55631,55503,55045,55056)),
    
  startle = list(
    transformation = function (x) log10(x),
    covariates       = c("sex","d50bw","ppibox1","ppibox2","ppibox3"),
    outliers         = c(55573,55043,55553,55292,55648,55866,55694,55693, 
                         55184,57029)))
