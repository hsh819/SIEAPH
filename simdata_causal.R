# Synthetic Dataset
# cite from Wei Sun, Pengyuan Wang, Dawei Yin, Jian Yang, and Yi Chang. Causal inference via sparse additive
# models with application to online advertising. In Proceedings of Twenty-Ninth AAAI Conference on
# Artificial Intelligence, pages 297–303, 2015.

# Define basis functions
f1 = function(x) -2 * sin(2 * x)
f2 = function(x) x^2 - 1/3
f3 = function(x) x - 0.5
f4 = function(x) exp(-x) - exp(-1) - 1
f5 = function(x) (x - 0.5)^2 + 2
f6 = function(x) as.numeric(x > 0)
f7 = function(x) exp(-x)
f8 = function(x) cos(x)
f9 = function(x) x^2
f10 = function(x) x

# Function to calculate T
cal_T = function(x) {
    T = f1(x[1]) + f2(x[2]) + f3(x[3]) + f4(x[4]) + f5(x[5])
    T = as.numeric(T > 0)
    return(T)
}

# Function to calculate Y
cal_Y = function(x) {
    f6(x[1]) + f7(x[2]) + f8(x[3]) + f9(x[4]) + f10(x[5])
}


simdata_synthetic = function(N = 1000, D = 100) {
  # Generate features independently
  X = sapply(1:D, function(i) rnorm(N, 0, 1))
  
  # Calculate treatment variable T
  T = apply(X[, 1:5], 1, cal_T)
  
  # Calculate outcomes Y
  Y_mean = apply(X[, 1:5], 1, cal_Y) + T 
  Y = Y_mean + rnorm(N, 0, 1)
   
  # Return the dataset
  dat = data.frame(X, T, Y)
  attr(dat, "Truth value of ATT") = 1
  return(dat)
}

# Generate the dataset
# set.seed(123) 
# dat = simdata_synthetic()



# Austin, Peter C. 2009. “Type I Error Rates, Coverage of Confidence Intervals, and Variance Estimation in Propensity-Score Matched Analyses.” 
# The International Journal of Biostatistics 5 (1). https://doi.org/10.2202/1557-4679.1146.
#Generating data similar to Austin (2009) for demonstrating treatment effect estimation

simdata_causal = function(N = 2000) {
  # Generated 9 covariates
  # x9 is independent of treatment and outcome
  X = sapply(1:9, function(i) rnorm(N, 0, 1))
  
  # 6 related to treatment selection: x1, x2, x4, x5, x7, x8
  LPT = log(2)*X[, 1] + log(1.5)*X[, 2] + log(2)*X[, 4] + log(1.5)*X[, 5] + log(2)*X[, 7] + log(1.5)*X[, 8] - 0.4
  PT = plogis(LPT)
  T = rbinom(N, 1, PT)  # approximately 40% of subjects would be exposed to the treatment

  # 6 related to outcome: x1, x2, x3, x4, x5, x6
  Y = 2*T + 2*X[, 1] + 2*X[, 2] + 2*X[, 3] + 1*X[, 4] + 1*X[, 5] + 1*X[, 6] + rnorm(N, 0, 5)
  
  # Return the dataset
  dat = data.frame(X, T, Y)
  attr(dat, "Truth value of ATT") = 2
  return(dat)  
  
}

# Generate the dataset
# set.seed(123)
# dat = simdata_causal()

			 