#' Fitting Weighted Quantile Sum regression models
#'
#' Fits Weighted Quantile Sum (WQS) regressions for continuous, binomial with random effect.
#'
#' @param formula An object of class \code{formula} specifying the relationship to be tested. The \code{wqs}
#' term must be included in \code{formula}, e.g. \code{y ~ wqs + ...}. To test for an interaction term with
#' a continuous variable \code{a} or for a quadratic term we can specify the \code{formula} as below:
#' \code{y ~ wqs*a + ...} and \code{y ~ wqs + I(wqs^2) + ...}, respectively.
#' @param data The \code{data.frame} containing the variables to be included in the model.
#' @param weights an optional vector of weights to be used in the fitting process.
#' Should be \code{NULsL} or a numeric vector.
#' @param mix_name A character vector listing the variables contributing to a mixture effect.
#' @param stratified The character name of the variable for which you want to stratify for.
#' It has to be a \code{factor}.
#' @param valid_var A character value containing the name of the variable that identifies the validation
#' and the training dataset. You previously need to create a variable in the dataset which is equal to 1
#' for the observations you want to include in the validation dataset, equal to 0 for the observation
#' you want to include in the training dataset (use 0 also for the validation dataset if you want to train and
#' validate the model on the same data) and equal to 2 if you want to keep part of the data for the
#' predictive model.
#' @param b Number of bootstrap samples used in parameter estimation.
#' @param b1_pos A logical value that determines whether weights are derived from models where the beta
#' values were positive or negative.
#' @param b1_constr A logial value that determines whether to apply positive (if \code{b1_pos = TRUE}) or
#' negative (if \code{b1_pos = FALSE}) constraints in the optimization function for the weight estimation.
#' @param q An \code{integer} to specify how mixture variables will be ranked, e.g. in quartiles
#' (\code{q = 4}), deciles (\code{q = 10}), or percentiles (\code{q = 100}). If \code{q = NULL} then
#' the values of the mixture variables are taken (these must be standardized).
#' @param validation Percentage of the dataset to be used to validate the model. If
#' \code{validation = 0} then the test dataset is used as validation dataset too.
#' @param family A character value that allows to decide for the glm: \code{gaussian} for linear regression,
#' \code{binomial} for logistic regression \code{"multinomial"} for multinomial regression,
#' \code{poisson} for Poisson regression, \code{quasipoisson} for quasi-Poisson regression,
#' \code{"negbin"} for negative binomial regression.
#' @param seed An \code{integer} value to fix the seed, if it is equal to \code{NULL} no seed is chosen.
#' @param pred Percentage of the dataset to be used for the predictive model. If \code{pred = 0} then no
#' predicitve model is going to be built.
#' @param plan_strategy A character value that allows to choose the evaluation strategies for the
#' \code{plan} function. You can choose among "sequential", "transparent", "multisession", "multicore",
#' "multiprocess", "cluster" and "remote" (see \code{\link[future]{plan}} help page for more details).
#' @param control The control list of optimization parameters. See \code{\link[Rsolnp]{solnp}} for details.
#' @param lp The lambda parameter that add a penlization term when we want to constrain in the negative direction.
#' This is an alternative to \code{b1_constr = TRUE}.
#' @param ln The lambda parameter that add a penlization term when we want to constrain in the positive direction.
#' This is an alternative to \code{b1_constr = TRUE}.
#'
#' @details
#' \code{gWQS} uses the \code{glm} function in the \bold{stats} package to fit the linear, logistic,
#' the Poisson and the quasi-Poisson regression, while the \code{glm.nb} function from the \bold{MASS}
#' package is used to fit the negative binomial regression respectively. The \code{nlm} function from
#' the \bold{stats} package was used to optimize the log-likelihood of the multinomial regression.\cr
#'
#' The \code{\link[Rsolnp]{solnp}} optimization function is used to estimate the weights at each
#' bootstrap step.\cr
#'
#' The \code{seed} argument specifies a fixed seed through the \code{\link[base]{set.seed}} function.\cr
#'
#' @return \code{gwqs} return the results of the WQS regression as well as many other objects and datasets.
#'
#' \item{fit}{The object that summarizes the output of the WQS model, reflecting a
#' linear, logistic, multinomial, Poisson, quasi-Poisson or negative binomial regression
#' depending on how the \code{family} parameter was specified.
#' The summary function can be used to call and print fit data (not for multinomial regression).}
#' \item{conv}{Indicates whether the solver has converged (0) or not (1 or 2).}
#' \item{bres}{Matrix of estimated weights, mixture effect parameter estimates and the associated
#' standard errors, statistics and p-values estimated for each bootstrap iteration.}
#' \item{wqs}{Vector containing the wqs index for each subject.}
#' \item{q_i}{List of the cutoffs used to divide in quantiles the variables in the mixture}
#' \item{bindex}{List of vectors containing the \code{rownames} of the subjects included in each
#' bootstrap dataset.}
#' \item{tindex}{Vector containing the rows used to estimate the weights in each bootstrap.}
#' \item{vindex}{Vector containing the rows used to estimate the parameters of the final model.}
#' \item{final_weights}{\code{data.frame} containing the final weights associated to each chemical.}
#' \item{y_wqs_df}{\code{data.frame} containing the dependent variable values adjusted for the
#' residuals of a fitted model adjusted for covariates (original values when \code{family = binomial}
#' or \code{"multinomial"}) and the wqs index estimated values.}
#' \item{df_pred}{\code{data.frame} containing the variables to print the ROC curve. It is generated only
#' when \code{pred > 0}}
#' \item{pindex}{Vector containing the subjects used for prediction. It is generated only when \code{pred > 0}}
#'
#' @author
#' Stefano Renzetti, Paul Curtin, Allan C Just, Ghalib Bello, Chris Gennings
#'
#' @references
#' Renzetti S, Gennings C, Curtin PC. 2019. gWQS: An R Package for Linear and Generalized Weighted
#' Quantile Sum (WQS) Regression. Journal of Statistical Software.\cr
#'
#' Carrico C, Gennings C, Wheeler D, Factor-Litvak P. Characterization of a weighted quantile sum
#' regression for highly correlated data in a risk analysis setting. J Biol Agricul Environ Stat.
#' 2014:1-21. ISSN: 1085-7117. DOI: 10.1007/ s13253-014-0180-3.
#' \url{http://dx.doi.org/10.1007/s13253-014-0180-3}.\cr
#'
#' Czarnota J, Gennings C, Colt JS, De Roos AJ, Cerhan JR, Severson RK, Hartge P, Ward MH,
#' Wheeler D. 2015. Analysis of environmental chemical mixtures and non-Hodgkin lymphoma risk in the
#' NCI-SEER NHL study. Environmental Health Perspectives, DOI:10.1289/ehp.1408630.\cr
#'
#' Czarnota J, Gennings C, Wheeler D. 2015. Assessment of weighted quantile sum regression for modeling
#' chemical mixtures and cancer risk. Cancer Informatics,
#' 2015:14(S2) 159-171 DOI: 10.4137/CIN.S17295.\cr
#'
#' Brunst KJ, Sanchez Guerra M, Gennings C, et al. Maternal Lifetime Stress and Prenatal Psychological
#' Functioning and Decreased Placental Mitochondrial DNA Copy Number in the PRISM Study.
#' Am J Epidemiol. 2017;186(11):1227-1236. doi:10.1093/aje/kwx183.\cr
#'
#' @examples
#' # we save the names of the mixture variables in the variable "toxic_chems"
#' toxic_chems = c("log_LBX074LA", "log_LBX099LA", "log_LBX105LA", "log_LBX118LA",
#' "log_LBX138LA", "log_LBX153LA", "log_LBX156LA", "log_LBX157LA", "log_LBX167LA",
#' "log_LBX170LA", "log_LBX180LA", "log_LBX187LA", "log_LBX189LA", "log_LBX194LA",
#' "log_LBX196LA", "log_LBX199LA", "log_LBXD01LA", "log_LBXD02LA", "log_LBXD03LA",
#' "log_LBXD04LA", "log_LBXD05LA", "log_LBXD07LA", "log_LBXF01LA", "log_LBXF02LA",
#' "log_LBXF03LA", "log_LBXF04LA", "log_LBXF05LA", "log_LBXF06LA", "log_LBXF07LA",
#' "log_LBXF08LA", "log_LBXF09LA", "log_LBXPCBLA", "log_LBXTCDLA", "log_LBXHXCLA")
#'
#' # To run a linear model and save the results in the variable "results". This linear model
#' # (family = gaussian) will rank/standardize variables in quartiles (q = 4), perform a
#' # 40/60 split of the data for training/validation (validation = 0.6), and estimate weights
#' # over 2 bootstrap samples (b = 2; in practical applications at least 100 bootstraps
#' # should be used). Weights will be derived from mixture effect parameters that are positive
#' # (b1_pos = TRUE). A unique seed was specified (seed = 2016) so this model will be
#' # reproducible, and plots describing the variable weights and linear relationship will be
#' # generated as output (plots = TRUE). In the end tables describing the weights values and
#' # the model parameters with the respectively statistics are generated in the plots window
#' # (tables = TRUE):
#' results = gwqs(y ~ wqs, mix_name = toxic_chems, data = wqs_data, q = 4, validation = 0.6,
#'                b = 2, b1_pos = TRUE, b1_constr = FALSE, family = gaussian, id='ID', seed = 2016
#'                )
#'
#' # to test the significance of the covariates
#' summary(results$fit)
#'

#' @import Rsolnp
#' @import stats
#'
#' @importFrom future plan future value
#' @importFrom future.apply future_lapply
#'
#' @export

# results = gwqs(glucose ~ wqs, mix_name = pollutants, data = lst$dat, q = 10, validation = 0.6,
                # b = 200, b1_pos = TRUE, b1_constr = FALSE, family = gaussian, id='ID', seed = 123, plan_strategy = "multisession")  


library(future)
library(future.apply)
library(Rsolnp)
library(lme4)
             

gwqs_re <- function(formula, data, weights, mix_name, stratified, valid_var, b = 100,
                 b1_pos = TRUE, b1_constr = FALSE, q = 4, validation = 0.6,
                 family = gaussian, id=NULL, seed = NULL, pred = 0, plan_strategy = "sequential",
                 control = list(rho = 1, outer.iter = 400, inner.iter = 800, delta = 1.0e-7, tol = 1e-8, trace = 0),
                 lp = 0, ln = 0)
{

	if(is.character(family)){
	  if(!family %in% c("gaussian", "binomial") ){
		stop("Only gaussian and binomial family supported")
	  } 
	  family <- get(family, mode = "function", envir = parent.frame())  
	}

	if(is.function(family)) {
	  tmp <- family()
	  if(!tmp$family %in% c("gaussian", "binomial")) {
		stop("Only gaussian and binomial family supported")
	  }
	  family <- tmp  
	}

	if(is.null(family$family)) {
	  stop("Family not recognized")  
	}

  if(!any(grepl("wqs", rownames(attr(terms(formula), "factors"))))) stop("'formula' must contain 'wqs' term: e.g. y ~ wqs + ...")
 
  allvars = all.vars(formula)
  data$wqs = rnorm(nrow(data), 0, 1)  # 初始化wqs
  data <- data[, c(allvars, mix_name, id)]
  dtf <- na.omit(data)
  
  # select rows of the original dataset based on rows selected with na.action

  # defining quantile variables
  if (is.null(q)) {Q = as.matrix(dtf[, mix_name]); q_i = NULL}
  else {
    q_f = quantile_f(dtf, mix_name, q)
    Q = q_f$Q
    q_i = q_f$q_i
  }

  # create stratified variables if stratified is not NULL
  m <- 0
  if (m){
    if(is.null(q)) stop("q must be different from NULL if you want to stratify for a categorical variable")
    strtfd_out = stratified_f(Q, dtf, stratified, mix_name)
    Q <- strtfd_out$Q
    mix_name = strtfd_out$mix_name
  }

  if(!is.null(seed) & !is.numeric(seed)) stop("seed must be numeric or NULL")
  if(!is.null(seed)) set.seed(seed)

  N <- nrow(dtf)

  # splitting the dataset
  m <- 0
  rindex = create_rindex(dtf, N, validation, pred, valid_var, m, family)

  # parameters estimation and model fitting
  m <- 0
  if(m[1]) weights <- dtf[, weights, drop = FALSE]
  else  dtf$weights <- weights <- rep(1, N)
  par_model <- par.modl.est(dtf[rindex$it,], Q[rindex$it,], formula, weights[rindex$it], b, b1_pos, b1_constr, family, id, plan_strategy, control, lp, ln)
  
  bres = par_model$bres
  conv = par_model$conv
  bindex = par_model$bindex
  nfuneval = par_model$n_funeval
  n_levels = par_model$n_levels
  strata_names = par_model$strata_names

  mean_weight = mean_weight_f(mix_name, bres, conv, b1_pos, family, n_levels, strata_names)

  # fit the final model with the estimated weights
  wqs_model = model.fit(mean_weight, dtf[rindex$iv,], Q[rindex$iv,], family, id, formula, weights[rindex$iv], b1_pos)  # id

  # prediction
  if(!is.null(rindex$ip)) df_pred = predict_f(Q[rindex$ip,], dtf[rindex$ip,], mean_weight, wqs_model$m_f, formula)
  else df_pred = NULL

  # Plots
  data_plot <- data.frame(mix_name, mean_weight)

  # data_plot = data_plot[order(data_plot$mean_weight, decreasing = TRUE),]
  wqs_index = as.numeric(unlist(wqs_model$wqs))

  # creating the list of elements to return
  results = list(wqs_model$m_f, conv, bres, wqs_index, q_i, bindex, rindex$it, rindex$iv,
                 data_plot, df_pred, rindex$ip)
  names(results) = c("fit", "conv", "bres", "wqs", "q_i", "bindex", "tindex", "vindex",
                     "final_weights", "df_pred", "pindex")  

  return(results)
}

