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

# The function below produces Figure 7 (b).

# Function No. 1

ProbabilityofExtremePovertyAnalysisBvsct <- function(InitialCapital = c(1.5, 2, 3, 4), PovertyLine = 1, file = "/Users/josemiguelflorescontro/Desktop"){
  # InitialCapital indicates the initial capital. Note that, this is a vector that contains the values for which one wants to perform the analysis. Moreover, these values should match the ones in the Mathematica notebook: ProbabilityofExtremePovertyAnalysisBvsct.nb
  # file is the location where the .tex file that will be used to generate the plot is going to be saved.
  
  setwd(file)
    MyColors <- pal_jco()(length(InitialCapital) + 1) 
    MyLines <-seq(from = 1, to = length(InitialCapital) + 1, by = 1) 
    B <- vector()
    my.expressions <- vector()
    # The values below were generated using Mathematica. In particular, they are outputs of the notebook: ProbabilityofExtremePovertyAnalysisBvsct.nb.nb
    CT <- c(0.8, 1.05, 1.3, 1.55, 1.8, 2.05, 2.3, 2.55, 2.8, 3.05, 3.3, 3.55, 
            3.8, 4.05, 4.3, 4.55, 4.8, 5.05, 5.3, 5.55, 5.8, 6.05, 6.3, 6.55, 
            6.8, 7.05, 7.3, 7.55, 7.8, 8.05, 8.3, 8.55, 8.8, 9.05, 9.3, 9.55, 9.8)
    for (i in 1:length(InitialCapital)){
      if (i==1){
        my.expressions <-c(my.expressions, paste0("$x$", " = ", InitialCapital[i]))
        B <- c(37.0518, 30.5781, 26.5523, 23.7973, 21.7878, 20.2538, 19.0422, 
               18.0595, 17.2453, 16.5589, 15.9718, 15.4635, 15.0187, 14.626, 
               14.2765, 13.9633, 13.6808, 13.4247, 13.1912, 12.9775, 12.781, 
               12.5998, 12.4319, 12.2761, 12.1309, 11.9953, 11.8683, 11.7492, 
               11.6371, 11.5315, 11.4318, 11.3375, 11.2482, 11.1634, 11.0828, 
               11.0061, 10.9331)
        # This is the name of the .tex file that will be saved in the location specified in the function.
        tikz('ProbabilityofExtremePovertyAnalysisBvsct.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
        par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
        plot(CT, B, type = "l", lwd = 2, lty = MyLines[2], col = MyColors[2], xaxs = "i", yaxs = "i", xlab = "$c_{\\scaleto{{T}}{4pt}}$", ylab = "$B$", ylim = c(0, 25), xlim = c(0, 8), cex.lab = 1, cex.axis = 1, xaxt = "n")
        axis(1, at = c(0, 2, 4, 6, 8), labels = FALSE)
        mtext(c("0", "2", "4", "6", "8"), side = 1, line = 1, at = c(0, 2, 4, 6, 8), col = c(rep("black", 5)))
      }else if (i == 2){
        my.expressions <-c(my.expressions, paste0("$x$", " = ", InitialCapital[i]))
        B <- c(34.0726, 28.334, 24.7565, 22.3027, 20.5092, 19.1374, 18.052, 17.1702, 
               16.4384, 15.8207, 15.2915, 14.8328, 14.431, 14.0758, 13.7593, 
               13.4754, 13.2192, 12.9866, 12.7744, 12.58, 12.4012, 12.236, 12.083, 
               11.9408, 11.8083, 11.6844, 11.5684, 11.4594, 11.3569, 11.2602, 
               11.1688, 11.0824, 11.0004, 10.9226, 10.8487, 10.7783, 10.7111)
        lines(CT, B, type = "l", lty = MyLines[i + 1], lwd = 2, col = MyColors[i + 1])
      } else if (i == 3){
        my.expressions <-c(my.expressions, paste0("$x$", " = ", InitialCapital[i]))
        B <- c(30.1415, 25.3593, 22.366, 20.3054, 18.7944, 17.6353, 16.7156, 
               15.9665, 15.3436, 14.8165, 14.3641, 13.9713, 13.6265, 13.3212, 
               13.0489, 12.8042, 12.5829, 12.3819, 12.1983, 12.0299, 11.8747, 
               11.7313, 11.5983, 11.4746, 11.3592, 11.2512, 11.1499, 11.0547, 
               10.9651, 10.8805, 10.8006, 10.7248, 10.653, 10.5847, 10.5198, 
               10.4579, 10.3989)
        lines(CT, B, type = "l", lty = MyLines[i + 1], lwd = 2, col = MyColors[i + 1])
      } else{
        my.expressions <-c(my.expressions, paste0("$x$", " = ", InitialCapital[i]))
        B <- c(27.504, 23.3465, 20.7363, 18.9346, 17.6103, 16.5922, 15.7829, 
               15.1227, 14.5726, 14.1066, 13.7061, 13.3579, 13.0519, 12.7807, 
               12.5384, 12.3206, 12.1235, 11.9442, 11.7803, 11.6298, 11.4911, 
               11.3628, 11.2438, 11.1329, 11.0295, 10.9326, 10.8417, 10.7563, 
               10.6757, 10.5997, 10.5278, 10.4596, 10.3949, 10.3335, 10.275, 
               10.2192, 10.166)
        lines(CT, B, type = "l", lty = MyLines[i + 1], lwd = 2, col = MyColors[i + 1])
      }
    }
    legend("topright", inset = 0.02, legend = my.expressions, lty = MyLines[2:length(MyLines)], lwd = 2, col = MyColors[2:length(MyLines)], cex = 0.8)
    dev.off()
  }
  
ProbabilityofExtremePovertyAnalysisBvsct(InitialCapital = c(1.5, 2, 3, 4), PovertyLine = 1, file = "/Users/josemiguelflorescontro/Desktop")
  
  
  
