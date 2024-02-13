#################################################################################################################
###################### Article: The Role of Direct Capital Cash Transfers Towards Poverty #######################
########################### and Extreme Poverty Alleviation - An Omega Risk Process #############################
########################### Authors: Flores-Contró, José M. and Arnold, Séverine ################################
#################################################################################################################

################################################################################################################
############################################ LOADING PACKAGES ##################################################
################################################################################################################

# We load the packages.

library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"))
library(ggsci)

################################################################################################################
########################################## 1st SET OF FUNCTIONS ################################################
################################################################################################################

# The functions below are in charge of simulating the Compound Poisson Process in the case where we 
# have jumps that follow a Beta Distribution with parameters alpha and 1. Both approaches should
# lead to very similar results (differences just due to randomness, number of simulations, etc.)

# Two main approaches were followed. These approaches are explained:

# Approach 1: In this approach we consider the fact that the waiting times
# are exponentially distributed with parameter lambda. Therefore, the Approach 1
# basically generates the waiting times of the events according to an exponential 
# distribution with parameter lambda. Then, the times when the events occur are added
# to the time grid of the Euler-Maruyama method. Moreover, a vector with the jumps is created.

# Function No.1

CompoundPoissonProcessBetaApproach1 <- function(T, lambda, alpha){
  # T represents the total time of the simulation (termination).
  # lambda represents the intensity parameter of the Poisson process (lambda > 0).
  # alpha represents the shape parameter of the Beta(alpha, 1) distribution (alpha > 0).
  t1 <- 0
  while (sum(t1) <= T){
    t1 <- c(t1, rexp(n = 1, rate = lambda))
  }
  t2 <- cumsum(t1[1:length(t1)-1])
  t3 <- sort(c(t2[2:length(t2)]))
  z <- rep(1,length(t3))
  for (i in match(t2[2:length(t2)],t3)){
    z[i] <- rbeta(n = 1, shape1 = alpha, shape2 = 1)
  }
  result <- list(t2, t3, z)
  return(result)}

# Example (Uncomment the below statement):
# CompoundPoissonProcessBetaApproach1(T = 20, lambda = 1, alpha = 1)

# Function No.2

# Approach 2: In this approach we consider the fact that there will be N
# events during the time interval [0, T] where N follows a Poisson Distribution
# with parameter lambda * T. In fact, this approach first generates the number
# of events that will occur in [0, T]. Then, by considering that these events could
# happen with the same probability at any time in [0,T] (using a uniform distribution)
# the function generates the times in [0, T] when the events occur. Finally, a vector with 
# the jumps is created.

CompoundPoissonProcessBetaApproach2 <- function(T, lambda, alpha){
  # T represents the total time of the simulation (termination).
  # lambda represents the intensity parameter of the Poisson process (lambda > 0).
  # alpha represents the shape parameter of the Beta(alpha, 1) distribution (alpha > 0).
  N <- rpois(n = 1, lambda * T)
  U <- runif (n = N, min = 0, max = T) 
  t1 <- sort(U)
  z <- rep(1,length(t1))
  for (i in match(U,t1)){
    z[i] <- rbeta(n = 1, shape1 = alpha, shape2 = 1)
  }
  result <- list(U, t1, z)
  return(result)}

# Example (Uncomment the below statement):
# CompoundPoissonProcessBetaApproach2(T = 20, lambda = 1, alpha = 1)

################################################################################################################
########################################## 2nd SET OF FUNCTIONS ################################################
################################################################################################################

# The function below is in charge of simulating the Capital Model described in Section 2. We use the Euler 
# Maruyama Method (see, for example, (Kloeden and Eckhard, 1995)) to simulate the stochastic process.

# Function No.1

EulerMaruyamaMethod <- function(T = 20, InitialCapital = 4, PovertyLine = 1, a = 0.1, b = 1.4, c = 0.4, cc = 1, B = 3, lambda = 1, alpha = 1, Simulations = 10){
    # T represents the total time of the simulation (termination).
    # InitialCapital indicates the Initial Capital of a Household (i.e., X(0) = Initial Capital).
    # PovertyLine represents the critical capital x* (or poverty line).
    # a is the rate of consumption (0 < a < 1).
    # b is the rate of income generation (0 < b).
    # c is the rate of investment or savings (0 < c < 1).
    # cc is the capital cash transfer rate (0 < cc).
    # B is the capital barrier level (B > PovertyLine).
    # lambda represents the intensity parameter of the Poisson process (lambda > 0).
    # alpha represents the shape parameter of the Beta(alpha, 1) distribution (alpha > 0).
    # Simulations is the number of simulations to be performed.
  r <- (1 - a) * b * c 
  for (k in 1:Simulations){ 
    ResultCompoundPoissonBeta <- CompoundPoissonProcessBetaApproach1(T, lambda, alpha) # Note that, here, we are particularly using the first approach from the 1st SET OF FUNCTIONS.
    x <- vector(,length(ResultCompoundPoissonBeta[[2]]))
    for (i in 1:length(ResultCompoundPoissonBeta[[2]])){
      if (i == 1){
        x[i] <- InitialCapital
      } else{
        if (x[i-1] >= PovertyLine && x[i-1] < B){
          tauB <- (1/(r - cc)) * log((B + ((cc * B) - (r * PovertyLine))/(r - cc))/(x[i-1] + ((cc * B) - (r * PovertyLine))/(r - cc)))
          time <- min(ResultCompoundPoissonBeta[[2]][i] - ResultCompoundPoissonBeta[[2]][i-1], tauB)
          if (time == tauB){
            x[i] <- ((B - PovertyLine) * exp(r * (ResultCompoundPoissonBeta[[1]][i] - ResultCompoundPoissonBeta[[1]][i-1] - time)) + PovertyLine) * ResultCompoundPoissonBeta[[3]][i-1]
          } else{
            x[i] <- ((x[i-1] + (cc * B - r * PovertyLine)/(r - cc)) * exp((r - cc) * (ResultCompoundPoissonBeta[[1]][i] - ResultCompoundPoissonBeta[[1]][i-1])) - (cc * B - r * PovertyLine)/(r - cc)) * ResultCompoundPoissonBeta[[3]][i-1]
          }
        } else if (x[i-1] >= B) {
            x[i] <- ((x[i-1] - PovertyLine) * exp(r * (ResultCompoundPoissonBeta[[1]][i] - ResultCompoundPoissonBeta[[1]][i-1])) + PovertyLine) * ResultCompoundPoissonBeta[[3]][i-1]
        } else if (x[i-1] < PovertyLine) {
          tauB1 <- (1/(r - cc)) * log((B + ((cc * B) - (r * PovertyLine))/(r - cc))/(x[i-1] + ((cc * B) - (r * PovertyLine))/(r - cc)))
          tauPovertyLine <- (-1/cc) * log((PovertyLine - B)/(x[i-1] - B))
          time1 <- min(ResultCompoundPoissonBeta[[2]][i] - ResultCompoundPoissonBeta[[2]][i-1], tauPovertyLine)
          if (time1 == tauPovertyLine){
            tauB2 <- (1/(r - cc)) * log((B + ((cc * B) - (r * PovertyLine))/(r - cc))/(PovertyLine + ((cc * B) - (r * PovertyLine))/(r - cc)))
            time2 <- min(ResultCompoundPoissonBeta[[2]][i] - ResultCompoundPoissonBeta[[2]][i-1] - tauPovertyLine, tauB2)
            if (time2 == tauB2){
              x[i] <- ((B - PovertyLine) * exp(r * (ResultCompoundPoissonBeta[[1]][i] - ResultCompoundPoissonBeta[[1]][i-1] - tauPovertyLine - time2)) + PovertyLine) * ResultCompoundPoissonBeta[[3]][i-1]
            } else {
              x[i] <- ((PovertyLine + (cc * B - r * PovertyLine)/(r - cc)) * exp((r - cc) * (ResultCompoundPoissonBeta[[1]][i] - ResultCompoundPoissonBeta[[1]][i-1] - tauPovertyLine)) - (cc * B - r * PovertyLine)/(r - cc)) * ResultCompoundPoissonBeta[[3]][i-1]
            } 
          } else {
              x[i] <- ((x[i-1] - B) * exp(-cc * (ResultCompoundPoissonBeta[[2]][i] - ResultCompoundPoissonBeta[[2]][i-1])) + B) * ResultCompoundPoissonBeta[[3]][i-1]
            }
          }
      }
    }
  }
  result <- list(ResultCompoundPoissonBeta[[1]], ResultCompoundPoissonBeta[[2]], ResultCompoundPoissonBeta[[3]], x)
  return(result)
}

# Example (Uncomment the below statement):
# EulerMaruyamaMethod(T = 20, InitialCapital = 0.5, PovertyLine = 1, a = 0.1, b = 0.4, c = 0.4, cc = 10, B = 2, lambda = 1, alpha = 3, Simulations = 10)

################################################################################################################
########################################## 3rd SET OF FUNCTIONS ################################################
################################################################################################################

# The functions below are in charge of implementing the Monte Carlo Simulation Methodology described in 
# Section 5.1. There is a function for each type of extreme poverty rate function (constant and exponential).

# Function No.1

ExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta <- function(T = 20, InitialCapital = 1.1, PovertyLine = 1, a = 0.1, b = 1.4, c = 0.4, cc = 0.504, B = 2, lambda = 1, alpha = 3, omega = 1, Simulations = 1){
  # T represents the total time of the simulation (termination).
  # InitialCapital indicates the Initial Capital of a Household (i.e., X(0) = Initial Capital).
  # PovertyLine represents the critical capital x* (or poverty line).
  # a is the rate of consumption (0 < a < 1).
  # b is the rate of income generation (0 < b).
  # c is the rate of investment or savings (0 < c < 1).
  # cc is the capital cash transfer rate (0 < cc).
  # B is the capital barrier level (B > PovertyLine).
  # lambda represents the intensity parameter of the Poisson process (lambda > 0).
  # alpha represents the shape parameter of the Beta(alpha, 1) distribution (alpha > 0).
  # omega represents the constant extreme poverty rate function described in Section 4.1.1.1 (omega > 0).
  # Simulations is the number of simulations to be performed.
  
  # set.seed(31) # Set up seed (if needed for reproducible results).
  r <- (1-a) * b * c 
  cat(sprintf("r: %s\n", r))
  Exponent <- 0
  CollectionOfExponents <- vector()
  psi <- vector()
  count <- 0
  for (k in 1:Simulations){
    Trajectory <- EulerMaruyamaMethod(T = T, InitialCapital = InitialCapital, PovertyLine = PovertyLine, a = a, b = b, c = c, cc = cc, B = B, lambda = lambda, alpha = alpha, Simulations = 1)
    phi <- vector()
    for (i in which(Trajectory[[4]] < PovertyLine)){
        phi <- c(phi, omega * min(Trajectory[[1]][i + 1] - Trajectory[[1]][i], - (1/cc * log((PovertyLine - B)/(Trajectory[[4]][i] - B)))))
    }
      phipertraj <- sum(phi)
      psi <- c(psi, 1 - exp(-phipertraj))
  }
  ProbabilityOfExtremePoverty <- mean(psi)
  CollectionOfExponentsStandardDeviaition <- (psi - ProbabilityOfExtremePoverty)^(2)
  StandardDeviation <- sqrt(1/(Simulations - 1) * sum(CollectionOfExponentsStandardDeviaition))
  LowerIntervalProbabilityOfExtremePoverty <- max(ProbabilityOfExtremePoverty - 2.81/sqrt(Simulations) * StandardDeviation, 0)
  UpperIntervalProbabilityOfExtremePoverty <- min(ProbabilityOfExtremePoverty + 2.81/sqrt(Simulations) * StandardDeviation, 1)
  cat(sprintf("Lower Interval Probability of Extreme Poverty (From Simulations): %s\n", LowerIntervalProbabilityOfExtremePoverty))
  cat(sprintf("Upper Interval Probability of Extreme Poverty (From Simulations): %s\n", UpperIntervalProbabilityOfExtremePoverty))
  cat(sprintf("Probability of Extreme Poverty (From Simulations): %s\n", ProbabilityOfExtremePoverty))
  result <- c(LowerIntervalProbabilityOfExtremePoverty, ProbabilityOfExtremePoverty, UpperIntervalProbabilityOfExtremePoverty)
  return(result)
}

# Example (Uncomment the below statement):
# ExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta <- ExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta(T = 1000, InitialCapital = 1, PovertyLine = 1, a = 0.1, b = 6, c = 0.4, cc = 0.1, B = 1.7, lambda = 1, alpha = 0.9, omega = 0.5, Simulations = 1000)

# Function No.2

ExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta <- function(T = 20, InitialCapital = 1.1, PovertyLine = 1, a = 0.1, b = 1.4, c = 0.4, cc = 0.504, B = 2, lambda = 1, alpha = 3, beta = 1, Simulations = 1){
  # T represents the total time of the simulation (termination).
  # InitialCapital indicates the Initial Capital of a Household (i.e., X(0) = Initial Capital).
  # PovertyLine represents the critical capital x* (or poverty line).
  # a is the rate of consumption (0 < a < 1).
  # b is the rate of income generation (0 < b).
  # c is the rate of investment or savings (0 < c < 1).
  # cc is the capital cash transfer rate (0 < cc).
  # B is the capital barrier level (B > PovertyLine).
  # lambda represents the intensity parameter of the Poisson process (lambda > 0).
  # alpha represents the shape parameter of the Beta(alpha, 1) distribution (alpha > 0).
  # beta represents the parameter of the extreme poverty rate function described in Section 4.1.1.2 (beta > 0).
  # Simulations is the number of simulations to be performed.
  
  # set.seed(31) # Set up seed (if needed for reproducible results).
  r <- (1-a) * b * c 
  cat(sprintf("r: %s\n", r))
  Exponent <- 0
  CollectionOfExponents <- vector()
  psi <- vector()
  count <- 0
  for (k in 1:Simulations){
    Trajectory <- EulerMaruyamaMethod(T = T, InitialCapital = InitialCapital, PovertyLine = PovertyLine, a = a, b = b, c = c, cc = cc, B = B, lambda = lambda, alpha = alpha, Simulations = 1)
    phi <- vector()
    for (i in which(Trajectory[[4]] < PovertyLine)){
      phi <- c(phi, beta/(B * cc) * (cc * (min(Trajectory[[1]][i + 1], Trajectory[[1]][i] - 1/cc * log((PovertyLine - B)/(Trajectory[[4]][i] - B))) - Trajectory[[1]][i]) + log(B + (Trajectory[[4]][i] - B) * exp(cc * (Trajectory[[1]][i] - min(Trajectory[[1]][i + 1], Trajectory[[1]][i] - 1/cc * log((PovertyLine - B)/(Trajectory[[4]][i] - B)))))) - log(B + (Trajectory[[4]][i] - B) * exp(cc * (Trajectory[[1]][i] - Trajectory[[1]][i])))))
    }
    phipertraj <- sum(phi)
    psi <- c(psi, 1 - exp(-phipertraj))
  }
  ProbabilityOfExtremePoverty <- mean(psi)
  CollectionOfExponentsStandardDeviaition <- (psi - ProbabilityOfExtremePoverty)^(2)
  StandardDeviation <- sqrt(1/(Simulations - 1) * sum(CollectionOfExponentsStandardDeviaition))
  LowerIntervalProbabilityOfExtremePoverty <- max(ProbabilityOfExtremePoverty - 2.81/sqrt(Simulations) * StandardDeviation, 0)
  UpperIntervalProbabilityOfExtremePoverty <- min(ProbabilityOfExtremePoverty + 2.81/sqrt(Simulations) * StandardDeviation, 1)
  cat(sprintf("Lower Interval Probability of Extreme Poverty (From Simulations): %s\n", LowerIntervalProbabilityOfExtremePoverty))
  cat(sprintf("Upper Interval Probability of Extreme Poverty (From Simulations): %s\n", UpperIntervalProbabilityOfExtremePoverty))
  cat(sprintf("Probability of Extreme Poverty (From Simulations): %s\n", ProbabilityOfExtremePoverty))
  result <- c(LowerIntervalProbabilityOfExtremePoverty, ProbabilityOfExtremePoverty, UpperIntervalProbabilityOfExtremePoverty)
  return(result)
}

# Example (Uncomment the below statement):
# ExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta <- ExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta(T = 1000, InitialCapital = 2.3, PovertyLine = 1, a = 0.1, b = 4, c = 0.4, cc = 1, B = 2, lambda = 1, alpha = 0.8, beta = 0.02, Simulations = 100)

################################################################################################################
########################################## 4rd SET OF FUNCTIONS ################################################
################################################################################################################

# The functions below are in charge of generating the plots shown in Section 5. Again, there is a function for 
# each type of extreme poverty rate function (constant and exponential). Moreover, there is a function for each
# type of analysis (capital cash transfer rate, capital barrier level and extreme poverty rate function).

# Function No.1

PlotCashTransferRateParameterAnalysisExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta <- function(T = 1000, InitialCapital = seq(0, 5, length = 30), PovertyLine = 1, a = 0.1, b = 4, c = 0.4, cc = c(0.25, 0.5, 0.75, 1), B = 2, lambda = 1, alpha = 0.8, omega = 0.02, Simulations = 10000, file = "/Users/josemiguelflorescontro/Desktop"){
  # T represents the total time of the simulation (termination).
  # InitialCapital indicates the Initial Capital of a Household (i.e., X(0) = Initial Capital).
  # PovertyLine represents the critical capital x* (or poverty line).
  # a is the rate of consumption (0 < a < 1).
  # b is the rate of income generation (0 < b).
  # c is the rate of investment or savings (0 < c < 1).
  # cc is the capital cash transfer rate (0 < cc). Note that, this is a vector that contains the values for which one wants to perform the analysis.
  # B is the capital barrier level (B > PovertyLine).
  # lambda represents the intensity parameter of the Poisson process (lambda > 0).
  # alpha represents the shape parameter of the Beta(alpha, 1) distribution (alpha > 0).
  # omega represents the constant extreme poverty rate function described in Section 4.1.1.1 (omega > 0).
  # Simulations is the number of simulations to be performed.
  # file is the location where the .tex file that will be used to generate the plot is going to be saved.
  
  par(mgp = c(3,1,0), mar = c(5,4,4,2) + 0.1)
  setwd(file)
  TrappingProbability <- vector()
  LowerProbabilityOfExtremePoverty <- vector()
  ProbabilityOfExtremePoverty <- vector()
  UpperProbabilityOfExtremePoverty <- vector()
  MyColors <- pal_jco()(length(cc) + 1)
  MyLines <-seq(from = 1, to = length(cc) + 1, by = 1)
  my.expressions <- vector()
  for (i in 1:length(cc)){
    for (k in 1:length(InitialCapital)){
      Calculation <- ExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta(T = T, InitialCapital = InitialCapital[k], PovertyLine = PovertyLine, a = a, b = b, c = c, cc = cc[i], B = B, lambda = lambda, alpha = alpha, omega = omega, Simulations = Simulations)
      LowerProbabilityOfExtremePoverty[k] <- Calculation[[1]]
      ProbabilityOfExtremePoverty[k] <- Calculation[[2]]
      UpperProbabilityOfExtremePoverty[k] <- Calculation[[3]]
    }
    if (i==1){
      # This is the name of the tex file that will be saved in the location specified in the function.
      tikz('PlotCashTransferRateParameterAnalysisExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
      my.expressions <- c(paste0("$c_{\\scaleto{T}{4pt}}$", " = ", cc[i]))
      par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
      plot(InitialCapital, ProbabilityOfExtremePoverty, type = "l", lty = MyLines[2], lwd = 2, col = MyColors[2], xaxs = "i", yaxs = "i", xlab = "Initial Capital", ylab = "$\\hat{\\psi}^{\\scaleto{{\\fontfamily{qcr}\\selectfont EP}}{2.5pt}}(x)_{n}$", ylim = c(0, 1), xlim = c(0, max(InitialCapital)), cex.lab = 1, cex.axis = 1, xaxt = "n")
      polygon(c(InitialCapital, rev(InitialCapital)), c(LowerProbabilityOfExtremePoverty, rev(UpperProbabilityOfExtremePoverty)), col = adjustcolor(MyColors[2], alpha.f = 0.3), border = FALSE)
      abline(v = PovertyLine, lty = 2, lwd = 0.5, col = "red")
      abline(v = B, lty = 2, lwd = 0.5, col = "#009900")
      axis(1, at = c(0, 1, 2, 3, 4, 5), labels = FALSE)
      mtext(c("0", "$x^{*} = 1$", "$B = 2$", "3", "4", "5"), side = 1, line = 1, at = c(0, 1, 2, 3, 4, 5), col = c("black", "red", "#009900", rep("black", 3)))
    }else{
      my.expressions <- c(my.expressions, paste0("$c_{\\scaleto{T}{4pt}}$", " = ", cc[i]) )
      lines(InitialCapital, ProbabilityOfExtremePoverty, type = "l", lty = MyLines[i + 1], lwd = 2, col = MyColors[i + 1])
      polygon(c(InitialCapital, rev(InitialCapital)), c(LowerProbabilityOfExtremePoverty, rev(UpperProbabilityOfExtremePoverty)), col = adjustcolor(MyColors[i + 1], alpha.f = 0.3), border = FALSE)
    }
  }
  legend("topright", inset = 0.02, legend = my.expressions, lty = MyLines[2:length(MyLines)], lwd = 2, col = MyColors[2:length(MyColors)], cex = 0.7)
  dev.off()
}

# Example (Uncomment the below statement):
# PlotCashTransferRateParameterAnalysisExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta(T = 1000, InitialCapital = seq(0, 5, length = 30), PovertyLine = 1, a = 0.1, b = 4, c = 0.4, cc = c(0.25, 0.5, 0.75, 1), B = 2, lambda = 1, alpha = 0.8, omega = 0.02, Simulations = 10000, file = "/Users/josemiguelflorescontro/Desktop")

# Function No.2

PlotCapitalBarrierLevelParameterAnalysisExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta <- function(T = 1000, InitialCapital = seq(0, 5, length = 30), PovertyLine = 1, a = 0.1, b = 4, c = 0.4, cc = 0.25, B = c(1.0001, 2, 3, 4), lambda = 1, alpha = 0.8, omega = 0.02, Simulations = 10000, file = "/Users/josemiguelflorescontro/Desktop"){
  # T represents the total time of the simulation (termination).
  # InitialCapital indicates the Initial Capital of a Household (i.e., X(0) = Initial Capital).
  # PovertyLine represents the critical capital x* (or poverty line).
  # a is the rate of consumption (0 < a < 1).
  # b is the rate of income generation (0 < b).
  # c is the rate of investment or savings (0 < c < 1).
  # cc is the capital cash transfer rate (0 < cc).
  # B is the capital barrier level (B > PovertyLine). Note that, this is a vector that contains the values for which one wants to perform the analysis.
  # lambda represents the intensity parameter of the Poisson process (lambda > 0).
  # alpha represents the shape parameter of the Beta(alpha, 1) distribution (alpha > 0).
  # omega represents the constant extreme poverty rate function described in Section 4.1.1.1 (omega > 0). 
  # Simulations is the number of simulations to be performed.
  # file is the location where the .tex file that will be used to generate the plot is going to be saved.
  
  par(mgp = c(3,1,0), mar = c(5,4,4,2) + 0.1)
  setwd(file)
  TrappingProbability <- vector()
  LowerProbabilityOfExtremePoverty <- vector()
  ProbabilityOfExtremePoverty <- vector()
  UpperProbabilityOfExtremePoverty <- vector()
  MyColors <- pal_jco()(length(B) + 1) 
  MyLines <-seq(from = 1, to = length(B) + 1, by = 1) 
  my.expressions <- vector()
  for (i in 1:length(B)){
    for (k in 1:length(InitialCapital)){
      Calculation <- ExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta(T = T, InitialCapital = InitialCapital[k], PovertyLine = PovertyLine, a = a, b = b, c = c, cc = cc, B = B[i], lambda = lambda, alpha = alpha, omega = omega, Simulations = Simulations)
      LowerProbabilityOfExtremePoverty[k] <- Calculation[[1]]
      ProbabilityOfExtremePoverty[k] <- Calculation[[2]]
      UpperProbabilityOfExtremePoverty[k] <- Calculation[[3]]
    }
    if (i==1){
      # This is the name of the tex file that will be saved in the location specified in the function.
      tikz('PlotCapitalBarrierLevelParameterAnalysisExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
      my.expressions <- c(paste0("$B$", " $\\rightarrow$ ", " $x^{*+}$ "))
      par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
      plot(InitialCapital, ProbabilityOfExtremePoverty, type = "l", lty = MyLines[2], lwd = 2, col = MyColors[2], xaxs = "i", yaxs = "i", xlab = "Initial Capital", ylab = "$\\hat{\\psi}^{\\scaleto{{\\fontfamily{qcr}\\selectfont EP}}{2.5pt}}(x)_{n}$", ylim = c(0, 1), xlim = c(0, max(InitialCapital)), cex.lab = 1, cex.axis = 1, xaxt = "n")
      polygon(c(InitialCapital, rev(InitialCapital)), c(LowerProbabilityOfExtremePoverty, rev(UpperProbabilityOfExtremePoverty)), col = adjustcolor(MyColors[2], alpha.f = 0.3), border = FALSE)
      abline(v = PovertyLine, lty = 2, lwd = 0.5, col = "red")
      axis(1, at = c(0, 1, 2, 3, 4, 5), labels = FALSE)
      mtext(c("0", "$x^{*} = 1$", "2", "3", "4", "5"), side = 1, line = 1, at = c(0, 1, 2, 3, 4, 5), col = c("black", "red", rep("black", 4)))
    }else{
      my.expressions <- c(my.expressions, paste0("$B$", " = ", B[i]) )
      lines(InitialCapital, ProbabilityOfExtremePoverty, type = "l", lty = MyLines[i + 1], lwd = 2, col = MyColors[i + 1])
      polygon(c(InitialCapital, rev(InitialCapital)), c(LowerProbabilityOfExtremePoverty, rev(UpperProbabilityOfExtremePoverty)), col = adjustcolor(MyColors[i + 1], alpha.f = 0.3), border = FALSE)
    }
  }
  legend("bottomright", inset = 0.02, legend = my.expressions, lty = MyLines[2:length(MyLines)], lwd = 2, col = MyColors[2:length(MyColors)], cex = 0.7)
  dev.off()
}

# Example (Uncomment the below statement):
# PlotCapitalBarrierLevelParameterAnalysisExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta(T = 1000, InitialCapital = seq(0, 5, length = 30), PovertyLine = 1, a = 0.1, b = 4, c = 0.4, cc = 0.25, B = c(1.0001, 2, 3, 4), lambda = 1, alpha = 0.8, omega = 0.02, Simulations = 10000, file = "/Users/josemiguelflorescontro/Desktop")

# Function No.3

PlotOmegaParameterAnalysisExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta <- function(T = 1000, InitialCapital = seq(0.001, 5, length = 30), PovertyLine = 1, a = 0.1, b = 4, c = 0.4, cc = 0.25, B = 2, lambda = 1, alpha = 0.8, omega = c(0.02, 0.05, 0.09), Simulations = 10000, file = "/Users/josemiguelflorescontro/Desktop"){
  # T represents the total time of the simulation (termination).
  # InitialCapital indicates the Initial Capital of a Household (i.e., X(0) = Initial Capital).
  # PovertyLine represents the critical capital x* (or poverty line).
  # a is the rate of consumption (0 < a < 1).
  # b is the rate of income generation (0 < b).
  # c is the rate of investment or savings (0 < c < 1).
  # cc is the capital cash transfer rate (0 < cc).
  # B is the capital barrier level (B > PovertyLine).
  # lambda represents the intensity parameter of the Poisson process (lambda > 0).
  # alpha represents the shape parameter of the Beta(alpha, 1) distribution (alpha > 0).
  # omega represents the constant extreme poverty rate function described in Section 4.1.1.1 (omega > 0). Note that, this is a vector that contains the values for which one wants to perform the analysis.
  # Simulations is the number of simulations to be performed.
  # file is the location where the .tex file that will be used to generate the plot is going to be saved.
  
  par(mgp = c(3,1,0), mar = c(5,4,4,2) + 0.1)
  setwd(file)
  TrappingProbability <- vector()
  LowerProbabilityOfExtremePoverty <- vector()
  ProbabilityOfExtremePoverty <- vector()
  UpperProbabilityOfExtremePoverty <- vector()
  MyColors <- pal_jco()(length(omega) + 1) 
  MyLines <-seq(from = 1, to = length(omega) + 1, by = 1) 
  my.expressions <- vector()
  for (i in 1:length(omega)){
    for (k in 1:length(InitialCapital)){
      Calculation <- ExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta(T = T, InitialCapital = InitialCapital[k], PovertyLine = PovertyLine, a = a, b = b, c = c, cc = cc, B = B, lambda = lambda, alpha = alpha, omega = omega[i], Simulations = Simulations)
      LowerProbabilityOfExtremePoverty[k] <- Calculation[[1]]
      ProbabilityOfExtremePoverty[k] <- Calculation[[2]]
      UpperProbabilityOfExtremePoverty[k] <- Calculation[[3]]
    }
    if (i==1){
      # This is the name of the tex file that will be saved in the location specified in the function.
      tikz('PlotOmegaParameterAnalysisExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
      my.expressions <- c("Trapping Probability $\\left(\\omega_{c}\\equiv \\infty\\right)$", paste0("$\\omega_{c}$", " = ", omega[i]))
      par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
      # The values for the trapping probability (the numbers in the below plot statement come from the Mathematica notebook: MonteCarloSimulations.nb)
      # In particular, they are the output of the "sol = Chop[TrappingProbability[#] & /@ MyList]" statement.
      plot(c(1., 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 
             1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21, 1.22, 
             1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 1.33, 
             1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.4, 1.41, 1.42, 1.43, 1.44, 
             1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 
             1.56, 1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.63, 1.64, 1.65, 1.66, 
             1.67, 1.68, 1.69, 1.7, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77, 
             1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.87, 1.88, 
             1.89, 1.9, 1.91, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98, 1.99, 2., 
             2.01, 2.02, 2.03, 2.04, 2.05, 2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 
             2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.2, 2.21, 2.22, 
             2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.31, 2.32, 2.33, 
             2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 
             2.45, 2.46, 2.47, 2.48, 2.49, 2.5, 2.51, 2.52, 2.53, 2.54, 2.55, 
             2.56, 2.57, 2.58, 2.59, 2.6, 2.61, 2.62, 2.63, 2.64, 2.65, 2.66, 
             2.67, 2.68, 2.69, 2.7, 2.71, 2.72, 2.73, 2.74, 2.75, 2.76, 2.77, 
             2.78, 2.79, 2.8, 2.81, 2.82, 2.83, 2.84, 2.85, 2.86, 2.87, 2.88, 
             2.89, 2.9, 2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 3., 
             3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07, 3.08, 3.09, 3.1, 3.11, 
             3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.2, 3.21, 3.22, 
             3.23, 3.24, 3.25, 3.26, 3.27, 3.28, 3.29, 3.3, 3.31, 3.32, 3.33, 
             3.34, 3.35, 3.36, 3.37, 3.38, 3.39, 3.4, 3.41, 3.42, 3.43, 3.44, 
             3.45, 3.46, 3.47, 3.48, 3.49, 3.5, 3.51, 3.52, 3.53, 3.54, 3.55, 
             3.56, 3.57, 3.58, 3.59, 3.6, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 
             3.67, 3.68, 3.69, 3.7, 3.71, 3.72, 3.73, 3.74, 3.75, 3.76, 3.77, 
             3.78, 3.79, 3.8, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 
             3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99, 4., 
             4.01, 4.02, 4.03, 4.04, 4.05, 4.06, 4.07, 4.08, 4.09, 4.1, 4.11, 
             4.12, 4.13, 4.14, 4.15, 4.16, 4.17, 4.18, 4.19, 4.2, 4.21, 4.22, 
             4.23, 4.24, 4.25, 4.26, 4.27, 4.28, 4.29, 4.3, 4.31, 4.32, 4.33, 
             4.34, 4.35, 4.36, 4.37, 4.38, 4.39, 4.4, 4.41, 4.42, 4.43, 4.44, 
             4.45, 4.46, 4.47, 4.48, 4.49, 4.5, 4.51, 4.52, 4.53, 4.54, 4.55, 
             4.56, 4.57, 4.58, 4.59, 4.6, 4.61, 4.62, 4.63, 4.64, 4.65, 4.66, 
             4.67, 4.68, 4.69, 4.7, 4.71, 4.72, 4.73, 4.74, 4.75, 4.76, 4.77, 
             4.78, 4.79, 4.8, 4.81, 4.82, 4.83, 4.84, 4.85, 4.86, 4.87, 4.88, 
             4.89, 4.9, 4.91, 4.92, 4.93, 4.94, 4.95, 4.96, 4.97, 4.98, 4.99, 5.), c(0.96736, 0.966064, 0.964788, 0.963531, 0.962292, 0.96107, 0.959865, 
                                                                                     0.958675, 0.957502, 0.956343, 0.955199, 0.954069, 0.952953, 0.95185, 
                                                                                     0.95076, 0.949683, 0.948618, 0.947564, 0.946523, 0.945492, 0.944473, 
                                                                                     0.943464, 0.942466, 0.941477, 0.940499, 0.939531, 0.938572, 0.937622, 
                                                                                     0.936682, 0.93575, 0.934827, 0.933913, 0.933007, 0.932109, 0.931219, 
                                                                                     0.930337, 0.929463, 0.928596, 0.927737, 0.926885, 0.92604, 0.925202, 
                                                                                     0.924371, 0.923547, 0.922729, 0.921918, 0.921113, 0.920315, 0.919523, 
                                                                                     0.918737, 0.917957, 0.917182, 0.916414, 0.915651, 0.914894, 0.914143, 
                                                                                     0.913397, 0.912656, 0.91192, 0.91119, 0.910465, 0.909745, 0.909029, 
                                                                                     0.908319, 0.907614, 0.906913, 0.906217, 0.905526, 0.904839, 0.904156, 
                                                                                     0.903478, 0.902805, 0.902135, 0.90147, 0.90081, 0.900153, 0.8995, 
                                                                                     0.898852, 0.898207, 0.897566, 0.896929, 0.896296, 0.895667, 0.895042, 
                                                                                     0.89442, 0.893802, 0.893187, 0.892576, 0.891969, 0.891365, 0.890764, 
                                                                                     0.890167, 0.889573, 0.888982, 0.888395, 0.887811, 0.88723, 0.886652, 
                                                                                     0.886078, 0.885506, 0.884937, 0.884372, 0.883811, 0.883254, 0.882701, 
                                                                                     0.882151, 0.881606, 0.881063, 0.880525, 0.87999, 0.879459, 0.878931, 
                                                                                     0.878406, 0.877885, 0.877368, 0.876853, 0.876342, 0.875834, 0.875329, 
                                                                                     0.874828, 0.874329, 0.873834, 0.873342, 0.872852, 0.872366, 0.871883, 
                                                                                     0.871402, 0.870924, 0.870449, 0.869977, 0.869508, 0.869042, 0.868578, 
                                                                                     0.868116, 0.867658, 0.867202, 0.866748, 0.866298, 0.865849, 0.865403, 
                                                                                     0.86496, 0.864519, 0.86408, 0.863644, 0.863211, 0.862779, 0.86235, 
                                                                                     0.861923, 0.861499, 0.861076, 0.860656, 0.860238, 0.859822, 0.859409, 
                                                                                     0.858997, 0.858588, 0.858181, 0.857775, 0.857372, 0.856971, 0.856572, 
                                                                                     0.856175, 0.85578, 0.855386, 0.854995, 0.854606, 0.854218, 0.853832, 
                                                                                     0.853449, 0.853067, 0.852687, 0.852308, 0.851932, 0.851557, 0.851184, 
                                                                                     0.850813, 0.850443, 0.850075, 0.849709, 0.849345, 0.848982, 0.848621, 
                                                                                     0.848261, 0.847903, 0.847547, 0.847192, 0.846839, 0.846488, 0.846138, 
                                                                                     0.845789, 0.845442, 0.845097, 0.844753, 0.84441, 0.84407, 0.84373, 
                                                                                     0.843392, 0.843055, 0.84272, 0.842386, 0.842054, 0.841723, 0.841393, 
                                                                                     0.841065, 0.840738, 0.840413, 0.840089, 0.839766, 0.839444, 0.839124, 
                                                                                     0.838805, 0.838487, 0.838171, 0.837856, 0.837542, 0.837229, 0.836918, 
                                                                                     0.836608, 0.836299, 0.835991, 0.835684, 0.835379, 0.835075, 0.834772, 
                                                                                     0.83447, 0.834169, 0.83387, 0.833571, 0.833274, 0.832978, 0.832683, 
                                                                                     0.832389, 0.832096, 0.831804, 0.831514, 0.831224, 0.830936, 0.830648, 
                                                                                     0.830362, 0.830076, 0.829792, 0.829509, 0.829227, 0.828945, 0.828665, 
                                                                                     0.828386, 0.828107, 0.82783, 0.827554, 0.827279, 0.827004, 0.826731, 
                                                                                     0.826458, 0.826187, 0.825916, 0.825647, 0.825378, 0.82511, 0.824843, 
                                                                                     0.824577, 0.824312, 0.824048, 0.823785, 0.823523, 0.823261, 0.823001, 
                                                                                     0.822741, 0.822482, 0.822224, 0.821967, 0.821711, 0.821455, 0.821201, 
                                                                                     0.820947, 0.820694, 0.820442, 0.82019, 0.81994, 0.81969, 0.819441, 
                                                                                     0.819193, 0.818946, 0.818699, 0.818454, 0.818209, 0.817964, 0.817721, 
                                                                                     0.817478, 0.817236, 0.816995, 0.816755, 0.816515, 0.816276, 0.816038, 
                                                                                     0.815801, 0.815564, 0.815328, 0.815093, 0.814858, 0.814624, 0.814391, 
                                                                                     0.814159, 0.813927, 0.813696, 0.813466, 0.813236, 0.813007, 0.812779, 
                                                                                     0.812551, 0.812324, 0.812098, 0.811872, 0.811647, 0.811423, 0.811199, 
                                                                                     0.810976, 0.810754, 0.810532, 0.810311, 0.81009, 0.809871, 0.809651, 
                                                                                     0.809433, 0.809215, 0.808998, 0.808781, 0.808565, 0.808349, 0.808134, 
                                                                                     0.80792, 0.807706, 0.807493, 0.807281, 0.807069, 0.806858, 0.806647, 
                                                                                     0.806437, 0.806227, 0.806018, 0.80581, 0.805602, 0.805394, 0.805188, 
                                                                                     0.804981, 0.804776, 0.804571, 0.804366, 0.804162, 0.803959, 0.803756, 
                                                                                     0.803553, 0.803352, 0.80315, 0.802949, 0.802749, 0.802549, 0.80235, 
                                                                                     0.802152, 0.801953, 0.801756, 0.801559, 0.801362, 0.801166, 0.80097, 
                                                                                     0.800775, 0.80058, 0.800386, 0.800193, 0.8, 0.799807, 0.799615, 
                                                                                     0.799423, 0.799232, 0.799041, 0.798851, 0.798661, 0.798472, 0.798283, 
                                                                                     0.798095, 0.797907, 0.79772, 0.797533, 0.797346, 0.79716, 0.796975, 
                                                                                     0.796789, 0.796605, 0.796421, 0.796237, 0.796053, 0.795871, 0.795688, 
                                                                                     0.795506, 0.795325, 0.795143, 0.794963, 0.794783, 0.794603, 0.794423, 
                                                                                     0.794244, 0.794066), type = "l", lty = MyLines[1], lwd = 2, col = MyColors[1], xaxs = "i", yaxs = "i", xlab = "Initial Capital", ylab = "$\\hat{\\psi}^{\\scaleto{{\\fontfamily{qcr}\\selectfont EP}}{2.5pt}}(x)_{n}$", ylim = c(0, 1), xlim = c(0, max(InitialCapital)), cex.lab = 1, cex.axis = 1, xaxt = "n")
      lines(InitialCapital, ProbabilityOfExtremePoverty, type = "l", lty = MyLines[2], lwd = 2, col = MyColors[2])
      polygon(c(InitialCapital, rev(InitialCapital)), c(LowerProbabilityOfExtremePoverty, rev(UpperProbabilityOfExtremePoverty)), col = adjustcolor(MyColors[2], alpha.f = 0.3), border = FALSE)
      abline(v = PovertyLine, lty = 2, lwd = 0.5, col = "red")
      abline(v = B, lty = 2, lwd = 0.5, col = "#009900")
      axis(1, at = c(0, 1, 2, 3, 4, 5), labels = FALSE)
      mtext(c("0", "$x^{*} = 1$", "$B = 2$", "3", "4", "5"), side = 1, line = 1, at = c(0, 1, 2, 3, 4, 5), col = c("black", "red", "#009900", rep("black", 3)))
    }else{
      my.expressions <- c(my.expressions, paste0("$\\omega_{c}$", " = ", omega[i]) )
      lines(InitialCapital, ProbabilityOfExtremePoverty, type = "l", lty = MyLines[i + 1], lwd = 2, col = MyColors[i + 1])
      polygon(c(InitialCapital, rev(InitialCapital)), c(LowerProbabilityOfExtremePoverty, rev(UpperProbabilityOfExtremePoverty)), col = adjustcolor(MyColors[i + 1], alpha.f = 0.3), border = FALSE)
    }
  }
  legend("bottomright", inset = 0.02, legend = my.expressions, lty = MyLines, lwd = 2, col = MyColors, cex = 0.7)
  dev.off()
}

# Example (Uncomment the below statement):
# PlotOmegaParameterAnalysisExtremePovertyProbabilityConstantExtremePovertyRateUninsuredSimulationBeta(T = 1000, InitialCapital = seq(0.001, 5, length = 30), PovertyLine = 1, a = 0.1, b = 4, c = 0.4, cc = 0.25, B = 2, lambda = 1, alpha = 0.8, omega = c(0.02, 0.05, 0.09), Simulations = 10000, file = "/Users/josemiguelflorescontro/Desktop")

# Function No.4

PlotCashTransferRateParameterAnalysisExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta <- function(T = 1000, InitialCapital = seq(0.001, 5, length = 30), PovertyLine = 1, a = 0.1, b = 4, c = 0.4, cc = c(0.25, 0.5, 0.75, 1), B = 2, lambda = 1, alpha = 0.8, beta = 0.02, Simulations = 10000, file = "/Users/josemiguelflorescontro/Desktop"){
  # T represents the total time of the simulation (termination).
  # InitialCapital indicates the Initial Capital of a Household (i.e., X(0) = Initial Capital).
  # PovertyLine represents the critical capital x* (or poverty line).
  # a is the rate of consumption (0 < a < 1).
  # b is the rate of income generation (0 < b).
  # c is the rate of investment or savings (0 < c < 1).
  # cc is the capital cash transfer rate (0 < cc). Note that, this is a vector that contains the values for which one wants to perform the analysis.
  # B is the capital barrier level (B > PovertyLine).
  # lambda represents the intensity parameter of the Poisson process (lambda > 0).
  # alpha represents the shape parameter of the Beta(alpha, 1) distribution (alpha > 0).
  # beta represents the parameter of the extreme poverty rate function described in Section 4.1.1.2 (beta > 0).
  # Simulations is the number of simulations to be performed.
  # file is the location where the .tex file that will be used to generate the plot is going to be saved.
  
  par(mgp = c(3,1,0), mar = c(5,4,4,2) + 0.1)
  setwd(file)
  TrappingProbability <- vector()
  LowerProbabilityOfExtremePoverty <- vector()
  ProbabilityOfExtremePoverty <- vector()
  UpperProbabilityOfExtremePoverty <- vector()
  MyColors <- pal_jco()(length(cc) + 1) 
  MyLines <-seq(from = 1, to = length(cc) + 1, by = 1) 
  my.expressions <- vector()
  for (i in 1:length(cc)){
    for (k in 1:length(InitialCapital)){
      Calculation <- ExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta(T = T, InitialCapital = InitialCapital[k], PovertyLine = PovertyLine, a = a, b = b, c = c, cc = cc[i], B = B, lambda = lambda, alpha = alpha, beta = beta, Simulations = Simulations)
      LowerProbabilityOfExtremePoverty[k] <- Calculation[[1]]
      ProbabilityOfExtremePoverty[k] <- Calculation[[2]]
      UpperProbabilityOfExtremePoverty[k] <- Calculation[[3]]
    }
    if (i==1){
      # This is the name of the tex file that will be saved in the location specified in the function.
      tikz('PlotCashTransferRateParameterAnalysisExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
      my.expressions <- c(paste0("$c_{\\scaleto{T}{4pt}}$", " = ", cc[i]))
      par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
      plot(InitialCapital, ProbabilityOfExtremePoverty, type = "l", lty = MyLines[2], lwd = 2, col = MyColors[2], xaxs = "i", yaxs = "i", xlab = "Initial Capital", ylab = "$\\hat{\\psi}^{\\scaleto{{\\fontfamily{qcr}\\selectfont EP}}{2.5pt}}(x)_{n}$", ylim = c(0, 1), xlim = c(0, max(InitialCapital)), cex.lab = 1, cex.axis = 1, xaxt = "n")
      polygon(c(InitialCapital, rev(InitialCapital)), c(LowerProbabilityOfExtremePoverty, rev(UpperProbabilityOfExtremePoverty)), col = adjustcolor(MyColors[2], alpha.f = 0.3), border = FALSE)
      abline(v = PovertyLine, lty = 2, lwd = 0.5, col = "red")
      abline(v = B, lty = 2, lwd = 0.5, col = "#009900")
      axis(1, at = c(0, 1, 2, 3, 4, 5), labels = FALSE)
      mtext(c("0", "$x^{*} = 1$", "$B = 2$", "3", "4", "5"), side = 1, line = 1, at = c(0, 1, 2, 3, 4, 5), col = c("black", "red", "#009900", rep("black", 3)))
    }else{
      my.expressions <- c(my.expressions, paste0("$c_{\\scaleto{T}{4pt}}$", " = ", cc[i]) )
      lines(InitialCapital, ProbabilityOfExtremePoverty, type = "l", lty = MyLines[i + 1], lwd = 2, col = MyColors[i + 1])
      polygon(c(InitialCapital, rev(InitialCapital)), c(LowerProbabilityOfExtremePoverty, rev(UpperProbabilityOfExtremePoverty)), col = adjustcolor(MyColors[i + 1], alpha.f = 0.3), border = FALSE)
    }
  }
  legend("bottomright", inset = 0.02, legend = my.expressions, lty = MyLines[2:length(MyLines)], lwd = 2, col = MyColors[2:length(MyColors)], cex = 0.7)
  dev.off()
}

# Example (Uncomment the below statement):
# PlotCashTransferRateParameterAnalysisExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta(T = 1000, InitialCapital = seq(0.001, 5, length = 30), PovertyLine = 1, a = 0.1, b = 4, c = 0.4, cc = c(0.25, 0.5, 0.75, 1), B = 2, lambda = 1, alpha = 0.8, beta = 0.02, Simulations = 10000, file = "/Users/josemiguelflorescontro/Desktop")

# Function No.5

PlotCapitalBarrierLevelParameterAnalysisExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta <- function(T = 1000, InitialCapital = seq(0.001, 5, length = 30), PovertyLine = 1, a = 0.1, b = 4, c = 0.4, cc = 0.25, B = c(1.0001, 2, 3, 4), lambda = 1, alpha = 0.8, beta = 0.02, Simulations = 10000, file = "/Users/josemiguelflorescontro/Desktop"){
  # T represents the total time of the simulation (termination).
  # InitialCapital indicates the Initial Capital of a Household (i.e., X(0) = Initial Capital).
  # PovertyLine represents the critical capital x* (or poverty line).
  # a is the rate of consumption (0 < a < 1).
  # b is the rate of income generation (0 < b).
  # c is the rate of investment or savings (0 < c < 1).
  # cc is the capital cash transfer rate (0 < cc). 
  # B is the capital barrier level (B > PovertyLine). Note that, this is a vector that contains the values for which one wants to perform the analysis.
  # lambda represents the intensity parameter of the Poisson process (lambda > 0).
  # alpha represents the shape parameter of the Beta(alpha, 1) distribution (alpha > 0).
  # beta represents the parameter of the extreme poverty rate function described in Section 4.1.1.2 (beta > 0).
  # Simulations is the number of simulations to be performed.
  # file is the location where the .tex file that will be used to generate the plot is going to be saved.
  
  par(mgp = c(3,1,0), mar = c(5,4,4,2) + 0.1)
  setwd(file)
  TrappingProbability <- vector()
  LowerProbabilityOfExtremePoverty <- vector()
  ProbabilityOfExtremePoverty <- vector()
  UpperProbabilityOfExtremePoverty <- vector()
  MyColors <- pal_jco()(length(B) + 1)
  MyLines <-seq(from = 1, to = length(B) + 1, by = 1)
  my.expressions <- vector()
  for (i in 1:length(B)){
    for (k in 1:length(InitialCapital)){
      Calculation <- ExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta(T = T, InitialCapital = InitialCapital[k], PovertyLine = PovertyLine, a = a, b = b, c = c, cc = cc, B = B[i], lambda = lambda, alpha = alpha, beta = beta, Simulations = Simulations)
      LowerProbabilityOfExtremePoverty[k] <- Calculation[[1]]
      ProbabilityOfExtremePoverty[k] <- Calculation[[2]]
      UpperProbabilityOfExtremePoverty[k] <- Calculation[[3]]
    }
    if (i==1){
      # This is the name of the tex file that will be saved in the location specified in the function.
      tikz('PlotCapitalBarrierLevelParameterAnalysisExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
      my.expressions <- c(my.expressions, paste0("$B$", " $\\rightarrow$ ", " $x^{*+}$ "))
      par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
      plot(InitialCapital, ProbabilityOfExtremePoverty, type = "l", lty = MyLines[2], lwd = 2, col = MyColors[2], xaxs = "i", yaxs = "i", xlab = "Initial Capital", ylab = "$\\hat{\\psi}^{\\scaleto{{\\fontfamily{qcr}\\selectfont EP}}{2.5pt}}(x)_{n}$", ylim = c(0, 1), xlim = c(0, max(InitialCapital)), cex.lab = 1, cex.axis = 1, xaxt = "n")
      polygon(c(InitialCapital, rev(InitialCapital)), c(LowerProbabilityOfExtremePoverty, rev(UpperProbabilityOfExtremePoverty)), col = adjustcolor(MyColors[2], alpha.f = 0.3), border = FALSE)
      abline(v = PovertyLine, lty = 2, lwd = 0.5, col = "red")
      axis(1, at = c(0, 1, 2, 3, 4, 5), labels = FALSE)
      mtext(c("0", "$x^{*} = 1$", "2", "3", "4", "5"), side = 1, line = 1, at = c(0, 1, 2, 3, 4, 5), col = c("black", "red", rep("black", 4)))
    }else{
      my.expressions <- c(my.expressions, paste0("$B$", " = ", B[i]) )
      lines(InitialCapital, ProbabilityOfExtremePoverty, type = "l", lty = MyLines[i + 1], lwd = 2, col = MyColors[i + 1])
      polygon(c(InitialCapital, rev(InitialCapital)), c(LowerProbabilityOfExtremePoverty, rev(UpperProbabilityOfExtremePoverty)), col = adjustcolor(MyColors[i + 1], alpha.f = 0.3), border = FALSE)
    }
  }
  legend("bottomright", inset = 0.02, legend = my.expressions, lty = MyLines[2:length(MyLines)], lwd = 2, col = MyColors[2:length(MyColors)], cex = 0.7)
  dev.off()
}

# Example (Uncomment the below statement):
# PlotCapitalBarrierLevelParameterAnalysisExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta(T = 1000, InitialCapital = seq(0.001, 5, length = 30), PovertyLine = 1, a = 0.1, b = 4, c = 0.4, cc = 0.25, B = c(1.0001, 2, 3, 4), lambda = 1, alpha = 0.8, beta = 0.02, Simulations = 10000, file = "/Users/josemiguelflorescontro/Desktop")

# Function No.6

PlotBetaParameterAnalysisExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta <- function(T = 1000, InitialCapital = seq(0.001, 5, length = 30), PovertyLine = 1, a = 0.1, b = 4, c = 0.4, cc = 0.25, B = 2, lambda = 1, alpha = 0.8, beta = c(0.02, 0.05, 0.09), Simulations = 10000, file = "/Users/josemiguelflorescontro/Desktop"){
  # T represents the total time of the simulation (termination).
  # InitialCapital indicates the Initial Capital of a Household (i.e., X(0) = Initial Capital).
  # PovertyLine represents the critical capital x* (or poverty line).
  # a is the rate of consumption (0 < a < 1).
  # b is the rate of income generation (0 < b).
  # c is the rate of investment or savings (0 < c < 1).
  # cc is the capital cash transfer rate (0 < cc).
  # B is the capital barrier level (B > PovertyLine).
  # lambda represents the intensity parameter of the Poisson process (lambda > 0).
  # alpha represents the shape parameter of the Beta(alpha, 1) distribution (alpha > 0).
  # beta represents the parameter of the extreme poverty rate function described in Section 4.1.1.2 (beta > 0). Note that, this is a vector that contains the values for which one wants to perform the analysis.
  # Simulations is the number of simulations to be performed.
  # file is the location where the .tex file that will be used to generate the plot is going to be saved.
  
  par(mgp = c(3,1,0), mar = c(5,4,4,2) + 0.1)
  setwd(file)
  TrappingProbability <- vector()
  LowerProbabilityOfExtremePoverty <- vector()
  ProbabilityOfExtremePoverty <- vector()
  UpperProbabilityOfExtremePoverty <- vector()
  MyColors <- pal_jco()(length(beta) + 1)
  MyLines <-seq(from = 1, to = length(beta) + 1, by = 1) 
  my.expressions <- vector()
  for (i in 1:length(beta)){
    for (k in 1:length(InitialCapital)){
      Calculation <- ExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta(T = T, InitialCapital = InitialCapital[k], PovertyLine = PovertyLine, a = a, b = b, c = c, cc = cc, B = B, lambda = lambda, alpha = alpha, beta = beta[i], Simulations = Simulations)
      LowerProbabilityOfExtremePoverty[k] <- Calculation[[1]]
      ProbabilityOfExtremePoverty[k] <- Calculation[[2]]
      UpperProbabilityOfExtremePoverty[k] <- Calculation[[3]]
    }
    if (i==1){
      # This is the name of the tex file that will be saved in the location specified in the function.
      tikz('PlotBetaParameterAnalysisExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
      my.expressions <- c("Trapping Probability $\\left(\\omega_{c}\\equiv \\infty\\right)$", paste0("$\\beta$", " = ", beta[i]))
      par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
      # The values for the trapping probability (the numbers in the below plot statement come from the Mathematica notebook: MonteCarloSimulations.nb)
      # In particular, they are the output of the "sol = Chop[TrappingProbability[#] & /@ MyList]" statement.
      plot(c(1., 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 
             1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21, 1.22, 
             1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 1.33, 
             1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.4, 1.41, 1.42, 1.43, 1.44, 
             1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 
             1.56, 1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.63, 1.64, 1.65, 1.66, 
             1.67, 1.68, 1.69, 1.7, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77, 
             1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.87, 1.88, 
             1.89, 1.9, 1.91, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98, 1.99, 2., 
             2.01, 2.02, 2.03, 2.04, 2.05, 2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 
             2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.2, 2.21, 2.22, 
             2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.31, 2.32, 2.33, 
             2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 
             2.45, 2.46, 2.47, 2.48, 2.49, 2.5, 2.51, 2.52, 2.53, 2.54, 2.55, 
             2.56, 2.57, 2.58, 2.59, 2.6, 2.61, 2.62, 2.63, 2.64, 2.65, 2.66, 
             2.67, 2.68, 2.69, 2.7, 2.71, 2.72, 2.73, 2.74, 2.75, 2.76, 2.77, 
             2.78, 2.79, 2.8, 2.81, 2.82, 2.83, 2.84, 2.85, 2.86, 2.87, 2.88, 
             2.89, 2.9, 2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 3., 
             3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07, 3.08, 3.09, 3.1, 3.11, 
             3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.2, 3.21, 3.22, 
             3.23, 3.24, 3.25, 3.26, 3.27, 3.28, 3.29, 3.3, 3.31, 3.32, 3.33, 
             3.34, 3.35, 3.36, 3.37, 3.38, 3.39, 3.4, 3.41, 3.42, 3.43, 3.44, 
             3.45, 3.46, 3.47, 3.48, 3.49, 3.5, 3.51, 3.52, 3.53, 3.54, 3.55, 
             3.56, 3.57, 3.58, 3.59, 3.6, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 
             3.67, 3.68, 3.69, 3.7, 3.71, 3.72, 3.73, 3.74, 3.75, 3.76, 3.77, 
             3.78, 3.79, 3.8, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 
             3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99, 4., 
             4.01, 4.02, 4.03, 4.04, 4.05, 4.06, 4.07, 4.08, 4.09, 4.1, 4.11, 
             4.12, 4.13, 4.14, 4.15, 4.16, 4.17, 4.18, 4.19, 4.2, 4.21, 4.22, 
             4.23, 4.24, 4.25, 4.26, 4.27, 4.28, 4.29, 4.3, 4.31, 4.32, 4.33, 
             4.34, 4.35, 4.36, 4.37, 4.38, 4.39, 4.4, 4.41, 4.42, 4.43, 4.44, 
             4.45, 4.46, 4.47, 4.48, 4.49, 4.5, 4.51, 4.52, 4.53, 4.54, 4.55, 
             4.56, 4.57, 4.58, 4.59, 4.6, 4.61, 4.62, 4.63, 4.64, 4.65, 4.66, 
             4.67, 4.68, 4.69, 4.7, 4.71, 4.72, 4.73, 4.74, 4.75, 4.76, 4.77, 
             4.78, 4.79, 4.8, 4.81, 4.82, 4.83, 4.84, 4.85, 4.86, 4.87, 4.88, 
             4.89, 4.9, 4.91, 4.92, 4.93, 4.94, 4.95, 4.96, 4.97, 4.98, 4.99, 5.), c(0.96736, 0.966064, 0.964788, 0.963531, 0.962292, 0.96107, 0.959865, 
                                                                                     0.958675, 0.957502, 0.956343, 0.955199, 0.954069, 0.952953, 0.95185, 
                                                                                     0.95076, 0.949683, 0.948618, 0.947564, 0.946523, 0.945492, 0.944473, 
                                                                                     0.943464, 0.942466, 0.941477, 0.940499, 0.939531, 0.938572, 0.937622, 
                                                                                     0.936682, 0.93575, 0.934827, 0.933913, 0.933007, 0.932109, 0.931219, 
                                                                                     0.930337, 0.929463, 0.928596, 0.927737, 0.926885, 0.92604, 0.925202, 
                                                                                     0.924371, 0.923547, 0.922729, 0.921918, 0.921113, 0.920315, 0.919523, 
                                                                                     0.918737, 0.917957, 0.917182, 0.916414, 0.915651, 0.914894, 0.914143, 
                                                                                     0.913397, 0.912656, 0.91192, 0.91119, 0.910465, 0.909745, 0.909029, 
                                                                                     0.908319, 0.907614, 0.906913, 0.906217, 0.905526, 0.904839, 0.904156, 
                                                                                     0.903478, 0.902805, 0.902135, 0.90147, 0.90081, 0.900153, 0.8995, 
                                                                                     0.898852, 0.898207, 0.897566, 0.896929, 0.896296, 0.895667, 0.895042, 
                                                                                     0.89442, 0.893802, 0.893187, 0.892576, 0.891969, 0.891365, 0.890764, 
                                                                                     0.890167, 0.889573, 0.888982, 0.888395, 0.887811, 0.88723, 0.886652, 
                                                                                     0.886078, 0.885506, 0.884937, 0.884372, 0.883811, 0.883254, 0.882701, 
                                                                                     0.882151, 0.881606, 0.881063, 0.880525, 0.87999, 0.879459, 0.878931, 
                                                                                     0.878406, 0.877885, 0.877368, 0.876853, 0.876342, 0.875834, 0.875329, 
                                                                                     0.874828, 0.874329, 0.873834, 0.873342, 0.872852, 0.872366, 0.871883, 
                                                                                     0.871402, 0.870924, 0.870449, 0.869977, 0.869508, 0.869042, 0.868578, 
                                                                                     0.868116, 0.867658, 0.867202, 0.866748, 0.866298, 0.865849, 0.865403, 
                                                                                     0.86496, 0.864519, 0.86408, 0.863644, 0.863211, 0.862779, 0.86235, 
                                                                                     0.861923, 0.861499, 0.861076, 0.860656, 0.860238, 0.859822, 0.859409, 
                                                                                     0.858997, 0.858588, 0.858181, 0.857775, 0.857372, 0.856971, 0.856572, 
                                                                                     0.856175, 0.85578, 0.855386, 0.854995, 0.854606, 0.854218, 0.853832, 
                                                                                     0.853449, 0.853067, 0.852687, 0.852308, 0.851932, 0.851557, 0.851184, 
                                                                                     0.850813, 0.850443, 0.850075, 0.849709, 0.849345, 0.848982, 0.848621, 
                                                                                     0.848261, 0.847903, 0.847547, 0.847192, 0.846839, 0.846488, 0.846138, 
                                                                                     0.845789, 0.845442, 0.845097, 0.844753, 0.84441, 0.84407, 0.84373, 
                                                                                     0.843392, 0.843055, 0.84272, 0.842386, 0.842054, 0.841723, 0.841393, 
                                                                                     0.841065, 0.840738, 0.840413, 0.840089, 0.839766, 0.839444, 0.839124, 
                                                                                     0.838805, 0.838487, 0.838171, 0.837856, 0.837542, 0.837229, 0.836918, 
                                                                                     0.836608, 0.836299, 0.835991, 0.835684, 0.835379, 0.835075, 0.834772, 
                                                                                     0.83447, 0.834169, 0.83387, 0.833571, 0.833274, 0.832978, 0.832683, 
                                                                                     0.832389, 0.832096, 0.831804, 0.831514, 0.831224, 0.830936, 0.830648, 
                                                                                     0.830362, 0.830076, 0.829792, 0.829509, 0.829227, 0.828945, 0.828665, 
                                                                                     0.828386, 0.828107, 0.82783, 0.827554, 0.827279, 0.827004, 0.826731, 
                                                                                     0.826458, 0.826187, 0.825916, 0.825647, 0.825378, 0.82511, 0.824843, 
                                                                                     0.824577, 0.824312, 0.824048, 0.823785, 0.823523, 0.823261, 0.823001, 
                                                                                     0.822741, 0.822482, 0.822224, 0.821967, 0.821711, 0.821455, 0.821201, 
                                                                                     0.820947, 0.820694, 0.820442, 0.82019, 0.81994, 0.81969, 0.819441, 
                                                                                     0.819193, 0.818946, 0.818699, 0.818454, 0.818209, 0.817964, 0.817721, 
                                                                                     0.817478, 0.817236, 0.816995, 0.816755, 0.816515, 0.816276, 0.816038, 
                                                                                     0.815801, 0.815564, 0.815328, 0.815093, 0.814858, 0.814624, 0.814391, 
                                                                                     0.814159, 0.813927, 0.813696, 0.813466, 0.813236, 0.813007, 0.812779, 
                                                                                     0.812551, 0.812324, 0.812098, 0.811872, 0.811647, 0.811423, 0.811199, 
                                                                                     0.810976, 0.810754, 0.810532, 0.810311, 0.81009, 0.809871, 0.809651, 
                                                                                     0.809433, 0.809215, 0.808998, 0.808781, 0.808565, 0.808349, 0.808134, 
                                                                                     0.80792, 0.807706, 0.807493, 0.807281, 0.807069, 0.806858, 0.806647, 
                                                                                     0.806437, 0.806227, 0.806018, 0.80581, 0.805602, 0.805394, 0.805188, 
                                                                                     0.804981, 0.804776, 0.804571, 0.804366, 0.804162, 0.803959, 0.803756, 
                                                                                     0.803553, 0.803352, 0.80315, 0.802949, 0.802749, 0.802549, 0.80235, 
                                                                                     0.802152, 0.801953, 0.801756, 0.801559, 0.801362, 0.801166, 0.80097, 
                                                                                     0.800775, 0.80058, 0.800386, 0.800193, 0.8, 0.799807, 0.799615, 
                                                                                     0.799423, 0.799232, 0.799041, 0.798851, 0.798661, 0.798472, 0.798283, 
                                                                                     0.798095, 0.797907, 0.79772, 0.797533, 0.797346, 0.79716, 0.796975, 
                                                                                     0.796789, 0.796605, 0.796421, 0.796237, 0.796053, 0.795871, 0.795688, 
                                                                                     0.795506, 0.795325, 0.795143, 0.794963, 0.794783, 0.794603, 0.794423, 
                                                                                     0.794244, 0.794066), type = "l", lty = MyLines[1], lwd = 2, col = MyColors[1], xaxs = "i", yaxs = "i", xlab = "Initial Capital", ylab = "$\\hat{\\psi}^{\\scaleto{{\\fontfamily{qcr}\\selectfont EP}}{2.5pt}}(x)_{n}$", ylim = c(0, 1), xlim = c(0, max(InitialCapital)), cex.lab = 1, cex.axis = 1, xaxt = "n")
      lines(InitialCapital, ProbabilityOfExtremePoverty, type = "l", lty = MyLines[2], lwd = 2, col = MyColors[2])
      polygon(c(InitialCapital, rev(InitialCapital)), c(LowerProbabilityOfExtremePoverty, rev(UpperProbabilityOfExtremePoverty)), col = adjustcolor(MyColors[2], alpha.f = 0.3), border = FALSE)
      abline(v = PovertyLine, lty = 2, lwd = 0.5, col = "red")
      abline(v = B, lty = 2, lwd = 0.5, col = "#009900")
      axis(1, at = c(0, 1, 2, 3, 4, 5), labels = FALSE)
      mtext(c("0", "$x^{*} = 1$", "$B = 2$", "3", "4", "5"), side = 1, line = 1, at = c(0, 1, 2, 3, 4, 5), col = c("black", "red", "#009900", rep("black", 3)))
    }else{
      my.expressions <- c(my.expressions, paste0("$\\beta$", " = ", beta[i]) )
      lines(InitialCapital, ProbabilityOfExtremePoverty, type = "l", lty = MyLines[i + 1], lwd = 2, col = MyColors[i + 1])
      polygon(c(InitialCapital, rev(InitialCapital)), c(LowerProbabilityOfExtremePoverty, rev(UpperProbabilityOfExtremePoverty)), col = adjustcolor(MyColors[i + 1], alpha.f = 0.3), border = FALSE)
    }
  }
  legend("bottomright", inset = 0.02, legend = my.expressions, lty = MyLines, lwd = 2, col = MyColors, cex = 0.7)
  dev.off()
}

# Example (Uncomment the below statement):
# PlotBetaParameterAnalysisExtremePovertyProbabilityExponentialExtremePovertyRateUninsuredSimulationBeta(T = 1000, InitialCapital = seq(0.001, 5, length = 30), PovertyLine = 1, a = 0.1, b = 4, c = 0.4, cc = 0.25, B = 2, lambda = 1, alpha = 0.8, beta = c(0.02, 0.05, 0.09), Simulations = 10000, file = "/Users/josemiguelflorescontro/Desktop")

################################################################################################################
############################################### REFERENCES #####################################################
################################################################################################################

# References: 
# 1 .- Kloeden, P. E. and P. Eckhard (1995). Numerical Solution of Stochastic Differential Equations. United 
# States of America: Springer-Verlag.

