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

# The function below produces Figure 7 (a).

# Function No. 1

TrappingProbabilityAnalysisBvsct <- function(InitialCapital = c(1.5, 2, 3, 4), file = "/Users/josemiguelflorescontro/Desktop"){
  # InitialCapital indicates the initial capital. Note that, this is a vector that contains the values for which one wants to perform the analysis. Moreover, these values should match the ones in the Mathematica notebook: TrappingProbabilityAnalysisBvsct.nb
  # file is the location where the .tex file that will be used to generate the plot is going to be saved.
  
    setwd(file)
    MyColors <- pal_jco()(length(InitialCapital) + 1) 
    MyLines <-seq(from = 1, to = length(InitialCapital) + 1, by = 1) 
    B <- vector()
    my.expressions <- vector()
    # The values below were generated using Mathematica. In particular, they are outputs of the notebook: TrappingProbabilityAnalysisBvsct.nb
    CT <- c(0.8, 1.05, 1.3, 1.55, 1.8, 2.05, 2.3, 2.55, 2.8, 3.05, 3.3, 3.55, 
            3.8, 4.05, 4.3, 4.55, 4.8, 5.05, 5.3, 5.55, 5.8, 6.05, 6.3, 6.55, 
            6.8, 7.05, 7.3, 7.55, 7.8, 8.05, 8.3, 8.55, 8.8, 9.05, 9.3, 9.55, 9.8)
    for (i in 1:length(InitialCapital)){
      if (i==1){
        my.expressions <-c(my.expressions, paste0("$x$", " = ", InitialCapital[i]))
        B <- c(421.778, 335.367, 282.051, 245.838, 219.616, 199.736, 184.136, 
               171.562, 161.206, 152.524, 145.139, 138.776, 133.237, 128.369, 
               124.056, 120.207, 116.75, 113.627, 110.793, 108.207, 105.838, 103.66, 
               101.65, 99.7893, 98.0613, 96.4521, 94.9497, 93.5436, 92.2248, 
               90.9852, 89.8178, 88.7162, 87.6751, 86.6894, 85.7547, 84.8672, 84.0233)
        # This is the name of the .tex file that will be saved in the location specified in the function.
        tikz('TrappingProbabilityAnalysisBvsct.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
        par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
        plot(CT, B, type = "l", lwd = 2, lty = MyLines[2], col = MyColors[2], xaxs = "i", yaxs = "i", xlab = "$c_{\\scaleto{{T}}{4pt}}$", ylab = "$B$", ylim = c(0, 300), xlim = c(0, 8), cex.lab = 1, cex.axis = 1, xaxt = "n")
        axis(1, at = c(0, 2, 4, 6, 8), labels = FALSE)
        mtext(c("0", "2", "4", "6", "8"), side = 1, line = 1, at = c(0, 2, 4, 6, 8), col = c(rep("black", 5)))
      }else if (i == 2){
        my.expressions <-c(my.expressions, paste0("$x$", " = ", InitialCapital[i]))
        B <- c(389.848, 311.293, 262.788, 229.819, 205.928, 187.802, 173.568, 
               162.087, 152.624, 144.687, 137.929, 132.105, 127.03, 122.568, 
               118.613, 115.081, 111.907, 109.039, 106.433, 104.056, 101.877, 
               99.8719, 98.0211, 96.3068, 94.7142, 93.2305, 91.8448, 90.5474, 
               89.3301, 88.1854, 87.1069, 86.089, 85.1266, 84.2151, 83.3506, 
               82.5294, 81.7482)
        lines(CT, B, type = "l", lty = MyLines[i + 1], lwd = 2, col = MyColors[i + 1])
      } else if (i == 3){
        my.expressions <-c(my.expressions, paste0("$x$", " = ", InitialCapital[i]))
        B <- c(348.428, 280.088, 237.833, 209.074, 188.205, 172.352, 159.888, 
               149.822, 141.516, 134.541, 128.596, 123.467, 118.994, 115.057, 
               111.563, 108.441, 105.633, 103.093, 100.784, 98.6752, 96.7411, 
               94.9604, 93.3153, 91.7905, 90.373, 89.0515, 87.8165, 86.6596, 
               85.5734, 84.5514, 83.588, 82.6782, 81.8176, 81.002, 80.2281, 79.4927, 
               78.7928)
        lines(CT, B, type = "l", lty = MyLines[i + 1], lwd = 2, col = MyColors[i + 1])
      } else{
        my.expressions <-c(my.expressions, paste0("$x$", " = ", InitialCapital[i]))
        B <- c(321.363, 259.709, 221.541, 195.532, 176.638, 162.268, 150.958, 
               141.815, 134.263, 127.914, 122.5, 117.823, 113.742, 110.146, 106.954, 
               104.098, 101.528, 99.2025, 97.0866, 95.153, 93.3785, 91.7438, 
               90.2327, 88.8314, 87.5279, 86.3123, 85.1756, 84.1102, 83.1095, 
               82.1676, 81.2793, 80.44, 79.6458, 78.8929, 78.1781, 77.4986, 76.8517)
        lines(CT, B, type = "l", lty = MyLines[i + 1], lwd = 2, col = MyColors[i + 1])
      }
    }
    legend("topright", inset = 0.02, legend = my.expressions, lty = MyLines[2:length(MyLines)], lwd = 2, col = MyColors[2:length(MyLines)], cex = 0.8)
    dev.off()
  }
  
TrappingProbabilityAnalysisBvsct(InitialCapital = c(1.5, 2, 3, 4), file = "/Users/josemiguelflorescontro/Desktop")
  
  
  
