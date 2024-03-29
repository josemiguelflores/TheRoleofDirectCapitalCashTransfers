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

SensitivityAnalysisforRateofConsumptionTrappingProbability <- function(InitialCapitals = c(1.3, 1.7, 4, 6), file = "/Users/josemiguelflorescontro/Desktop"){
  # InitialCapitals indicates the Initial Capital of a Household (i.e., X(0) = Initial Capital). Note that, this is a vector that contains the values for which one wants to perform the analysis. Moreover, these values should match the ones in the Mathematica notebook: EffectsofUnderlyingFactorsontheTrappingProbability.nb.
  # file is the location where the .tex file that will be used to generate the plot is going to be saved.
  
  setwd(file)
  MyColors <- pal_jco()(length(2)) 
  MyLines <- seq(from = 1, to = length(InitialCapitals), by = 1)
  TrappingProbabilityCashTransfer <- vector()
  TrappingProbabilityOriginal <- list()
  a <- vector()
  b <- vector()
  c <- vector()
  # The values below were generated using Mathematica. In particular, they are outputs of the notebook: EffectsofUnderlyingFactorsontheTrappingProbability.nb
  RateofConsumption <- c(0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 
                         0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 
                         0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 
                         0.031, 0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 0.04, 
                         0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.05, 
                         0.051, 0.052, 0.053, 0.054, 0.055, 0.056, 0.057, 0.058, 0.059, 0.06, 
                         0.061, 0.062, 0.063, 0.064, 0.065, 0.066, 0.067, 0.068, 0.069, 0.07, 
                         0.071, 0.072, 0.073, 0.074, 0.075, 0.076, 0.077, 0.078, 0.079, 0.08, 
                         0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, 
                         0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099, 0.1, 
                         0.101, 0.102, 0.103, 0.104, 0.105, 0.106, 0.107, 0.108, 0.109, 0.11, 
                         0.111, 0.112, 0.113, 0.114, 0.115, 0.116, 0.117, 0.118, 0.119, 0.12, 
                         0.121, 0.122, 0.123, 0.124, 0.125, 0.126, 0.127, 0.128, 0.129, 0.13, 
                         0.131, 0.132, 0.133, 0.134, 0.135, 0.136, 0.137, 0.138, 0.139, 0.14, 
                         0.141, 0.142, 0.143, 0.144, 0.145, 0.146, 0.147, 0.148, 0.149, 0.15, 
                         0.151, 0.152, 0.153, 0.154, 0.155, 0.156, 0.157, 0.158, 0.159, 0.16, 
                         0.161, 0.162, 0.163, 0.164, 0.165, 0.166, 0.167, 0.168, 0.169, 0.17, 
                         0.171, 0.172, 0.173, 0.174, 0.175, 0.176, 0.177, 0.178, 0.179, 0.18, 
                         0.181, 0.182, 0.183, 0.184, 0.185, 0.186, 0.187, 0.188, 0.189, 0.19, 
                         0.191, 0.192, 0.193, 0.194, 0.195, 0.196, 0.197, 0.198, 0.199, 0.2, 
                         0.201, 0.202, 0.203, 0.204, 0.205, 0.206, 0.207, 0.208, 0.209, 0.21, 
                         0.211, 0.212, 0.213, 0.214, 0.215, 0.216, 0.217)
  TrappingProbabilityCashTransfer <- c(0.881349, 0.88188, 0.882411, 0.882943, 0.883475, 0.884008, 0.884541, 
                                       0.885074, 0.885607, 0.88614, 0.886674, 0.887208, 0.887743, 0.888277, 
                                       0.888812, 0.889348, 0.889883, 0.890419, 0.890955, 0.891491, 0.892028, 
                                       0.892565, 0.893102, 0.893639, 0.894177, 0.894715, 0.895253, 0.895791, 
                                       0.89633, 0.896869, 0.897408, 0.897948, 0.898487, 0.899027, 0.899567, 
                                       0.900108, 0.900649, 0.901189, 0.901731, 0.902272, 0.902813, 0.903355, 
                                       0.903897, 0.90444, 0.904982, 0.905525, 0.906068, 0.906611, 0.907154, 
                                       0.907698, 0.908242, 0.908786, 0.90933, 0.909874, 0.910419, 0.910964, 
                                       0.911509, 0.912054, 0.9126, 0.913145, 0.913691, 0.914237, 0.914783, 
                                       0.91533, 0.915876, 0.916423, 0.91697, 0.917517, 0.918064, 0.918612, 
                                       0.919159, 0.919707, 0.920255, 0.920803, 0.921351, 0.9219, 0.922448, 
                                       0.922997, 0.923546, 0.924095, 0.924644, 0.925193, 0.925743, 0.926292, 
                                       0.926842, 0.927392, 0.927941, 0.928492, 0.929042, 0.929592, 0.930142, 
                                       0.930693, 0.931243, 0.931794, 0.932345, 0.932896, 0.933447, 0.933998, 
                                       0.934549, 0.9351, 0.935652, 0.936203, 0.936755, 0.937306, 0.937858, 
                                       0.938409, 0.938961, 0.939513, 0.940065, 0.940617, 0.941169, 0.941721, 
                                       0.942273, 0.942825, 0.943377, 0.943929, 0.944481, 0.945034, 0.945586, 
                                       0.946138, 0.94669, 0.947243, 0.947795, 0.948347, 0.9489, 0.949452, 
                                       0.950004, 0.950556, 0.951109, 0.951661, 0.952213, 0.952765, 0.953317, 
                                       0.953869, 0.954421, 0.954973, 0.955525, 0.956077, 0.956629, 0.957181, 
                                       0.957732, 0.958284, 0.958836, 0.959387, 0.959938, 0.96049, 0.961041, 
                                       0.961592, 0.962143, 0.962694, 0.963244, 0.963795, 0.964346, 0.964896, 
                                       0.965446, 0.965996, 0.966546, 0.967096, 0.967646, 0.968195, 0.968745, 
                                       0.969294, 0.969843, 0.970392, 0.97094, 0.971489, 0.972037, 0.972585, 
                                       0.973133, 0.973681, 0.974228, 0.974775, 0.975322, 0.975869, 0.976415, 
                                       0.976962, 0.977508, 0.978054, 0.978599, 0.979144, 0.979689, 0.980234, 
                                       0.980778, 0.981322, 0.981866, 0.98241, 0.982953, 0.983496, 0.984038, 
                                       0.984581, 0.985122, 0.985664, 0.986205, 0.986746, 0.987287, 0.987827, 
                                       0.988366, 0.988906, 0.989445, 0.989983, 0.990522, 0.991059, 0.991597, 
                                       0.992134, 0.99267, 0.993206, 0.993742, 0.994277, 0.994812, 0.995346, 
                                       0.99588, 0.996413, 0.996946, 0.997478, 0.99801, 0.998542, 0.999072)
  # This is the name of the .tex file that will be saved in the location specified in the function.
  tikz('SensitivityAnalysisforRateofConsumptionTrappingProbability.tex', standAlone = TRUE, width = 4, height = 4, packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}", "\\usepackage{amssymb}", "\\usepackage{scalerel}"))
  par(mgp = c(2.5, 1, 0), mar = c(3.5, 3.5, 1, 1) + 0.1)
  plot(RateofConsumption, TrappingProbabilityCashTransfer, type = "l", lwd = 2, lty = MyLines[1],  col = "blue", xaxs = "i", yaxs = "i", xlab = "$a$", ylab = "Trapping Probability", ylim = c(0, 1), xlim = c(round(0, digits = 0), 0.2), cex.lab = 1, cex.axis = 1)
  TrappingProbabilityOriginal <- c(0.893958, 0.894456, 0.894954, 0.895451, 0.895949, 0.896447, 0.896945, 
                                   0.897443, 0.897941, 0.89844, 0.898938, 0.899437, 0.899935, 0.900434, 
                                   0.900933, 0.901432, 0.901931, 0.90243, 0.902929, 0.903428, 0.903928, 
                                   0.904427, 0.904926, 0.905426, 0.905925, 0.906425, 0.906925, 0.907424, 
                                   0.907924, 0.908424, 0.908924, 0.909424, 0.909924, 0.910424, 0.910924, 
                                   0.911424, 0.911925, 0.912425, 0.912925, 0.913425, 0.913926, 0.914426, 
                                   0.914926, 0.915427, 0.915927, 0.916428, 0.916928, 0.917429, 0.917929, 
                                   0.918429, 0.91893, 0.91943, 0.919931, 0.920431, 0.920932, 0.921432, 
                                   0.921933, 0.922433, 0.922933, 0.923434, 0.923934, 0.924434, 0.924935, 
                                   0.925435, 0.925935, 0.926435, 0.926935, 0.927435, 0.927935, 0.928435, 
                                   0.928935, 0.929435, 0.929935, 0.930434, 0.930934, 0.931434, 0.931933, 
                                   0.932433, 0.932932, 0.933431, 0.93393, 0.934429, 0.934928, 0.935427, 
                                   0.935926, 0.936425, 0.936923, 0.937422, 0.93792, 0.938418, 0.938916, 
                                   0.939414, 0.939912, 0.94041, 0.940907, 0.941405, 0.941902, 0.942399, 
                                   0.942896, 0.943393, 0.94389, 0.944386, 0.944883, 0.945379, 0.945875, 
                                   0.946371, 0.946866, 0.947362, 0.947857, 0.948352, 0.948847, 0.949342, 
                                   0.949836, 0.95033, 0.950824, 0.951318, 0.951812, 0.952305, 0.952798, 
                                   0.953291, 0.953784, 0.954276, 0.954768, 0.95526, 0.955752, 0.956243, 
                                   0.956735, 0.957225, 0.957716, 0.958206, 0.958696, 0.959186, 0.959676, 
                                   0.960165, 0.960654, 0.961142, 0.961631, 0.962119, 0.962606, 0.963094, 
                                   0.96358, 0.964067, 0.964553, 0.965039, 0.965525, 0.96601, 0.966495, 
                                   0.96698, 0.967464, 0.967948, 0.968431, 0.968914, 0.969397, 0.969879, 
                                   0.970361, 0.970842, 0.971323, 0.971804, 0.972284, 0.972764, 0.973243, 
                                   0.973722, 0.974201, 0.974679, 0.975156, 0.975634, 0.97611, 0.976586, 
                                   0.977062, 0.977537, 0.978012, 0.978486, 0.97896, 0.979434, 0.979906, 
                                   0.980379, 0.98085, 0.981321, 0.981792, 0.982262, 0.982732, 0.983201, 
                                   0.983669, 0.984137, 0.984605, 0.985071, 0.985538, 0.986003, 0.986468, 
                                   0.986933, 0.987397, 0.98786, 0.988322, 0.988784, 0.989246, 0.989706, 
                                   0.990166, 0.990626, 0.991085, 0.991543, 0.992, 0.992457, 0.992913, 
                                   0.993368, 0.993823, 0.994277, 0.99473, 0.995183, 0.995634, 0.996085, 
                                   0.996536, 0.996985, 0.997434, 0.997882, 0.99833, 0.998776, 0.999222)
  lines(RateofConsumption, TrappingProbabilityOriginal, type = "l", lty = MyLines[1], lwd = 2, col = "red")
  for (k in 1:length(InitialCapitals)){
    if (k == 2){
      TrappingProbabilityCashTransfer <- c(0.823374, 0.824098, 0.824822, 0.825548, 0.826274, 0.827001, 0.827729, 
                                           0.828458, 0.829188, 0.829919, 0.830651, 0.831383, 0.832117, 0.832851, 
                                           0.833586, 0.834322, 0.835059, 0.835797, 0.836536, 0.837276, 0.838016, 
                                           0.838758, 0.8395, 0.840244, 0.840988, 0.841733, 0.842479, 0.843226, 
                                           0.843973, 0.844722, 0.845471, 0.846222, 0.846973, 0.847725, 0.848478, 
                                           0.849232, 0.849987, 0.850743, 0.851499, 0.852257, 0.853015, 0.853774, 
                                           0.854535, 0.855296, 0.856058, 0.85682, 0.857584, 0.858348, 0.859114, 
                                           0.85988, 0.860647, 0.861415, 0.862184, 0.862954, 0.863725, 0.864496, 
                                           0.865269, 0.866042, 0.866816, 0.867591, 0.868367, 0.869144, 0.869921, 
                                           0.8707, 0.871479, 0.872259, 0.87304, 0.873822, 0.874605, 0.875389, 
                                           0.876173, 0.876959, 0.877745, 0.878532, 0.87932, 0.880108, 0.880898, 
                                           0.881689, 0.88248, 0.883272, 0.884065, 0.884859, 0.885653, 0.886449, 
                                           0.887245, 0.888042, 0.88884, 0.889639, 0.890439, 0.891239, 0.89204, 
                                           0.892843, 0.893645, 0.894449, 0.895254, 0.896059, 0.896865, 0.897672, 
                                           0.89848, 0.899289, 0.900098, 0.900908, 0.901719, 0.902531, 0.903344, 
                                           0.904157, 0.904971, 0.905786, 0.906602, 0.907418, 0.908236, 0.909054, 
                                           0.909873, 0.910692, 0.911513, 0.912334, 0.913156, 0.913978, 0.914802, 
                                           0.915626, 0.916451, 0.917276, 0.918103, 0.91893, 0.919758, 0.920586, 
                                           0.921416, 0.922246, 0.923077, 0.923908, 0.92474, 0.925573, 0.926407, 
                                           0.927241, 0.928076, 0.928912, 0.929749, 0.930586, 0.931424, 0.932262, 
                                           0.933101, 0.933941, 0.934782, 0.935623, 0.936465, 0.937307, 0.938151, 
                                           0.938994, 0.939839, 0.940684, 0.94153, 0.942376, 0.943223, 0.944071, 
                                           0.944919, 0.945768, 0.946617, 0.947467, 0.948318, 0.949169, 0.950021, 
                                           0.950874, 0.951727, 0.95258, 0.953434, 0.954289, 0.955144, 0.956, 
                                           0.956856, 0.957713, 0.958571, 0.959429, 0.960287, 0.961146, 0.962005, 
                                           0.962865, 0.963726, 0.964587, 0.965448, 0.96631, 0.967173, 0.968036, 
                                           0.968899, 0.969763, 0.970627, 0.971491, 0.972357, 0.973222, 0.974088, 
                                           0.974954, 0.975821, 0.976688, 0.977556, 0.978423, 0.979292, 0.98016, 
                                           0.981029, 0.981898, 0.982768, 0.983638, 0.984508, 0.985379, 0.98625, 
                                           0.987121, 0.987993, 0.988864, 0.989736, 0.990609, 0.991481, 0.992354, 
                                           0.993227, 0.9941, 0.994974, 0.995847, 0.996721, 0.997595, 0.99847)
      lines(RateofConsumption, TrappingProbabilityCashTransfer, type = "l", lty = MyLines[k], lwd = 2, col = "blue")
      TrappingProbabilityOriginal <- c(0.835154, 0.83585, 0.836546, 0.837242, 0.83794, 0.838638, 0.839337, 
                                       0.840037, 0.840737, 0.841438, 0.84214, 0.842842, 0.843545, 0.844249, 
                                       0.844954, 0.845659, 0.846365, 0.847071, 0.847779, 0.848487, 0.849195, 
                                       0.849905, 0.850615, 0.851326, 0.852037, 0.852749, 0.853462, 0.854176, 
                                       0.85489, 0.855605, 0.856321, 0.857037, 0.857754, 0.858472, 0.85919, 
                                       0.859909, 0.860629, 0.861349, 0.86207, 0.862792, 0.863515, 0.864238, 
                                       0.864961, 0.865686, 0.866411, 0.867137, 0.867863, 0.86859, 0.869318, 
                                       0.870046, 0.870775, 0.871505, 0.872235, 0.872966, 0.873698, 0.87443, 
                                       0.875163, 0.875897, 0.876631, 0.877366, 0.878102, 0.878838, 0.879575, 
                                       0.880312, 0.88105, 0.881789, 0.882528, 0.883268, 0.884009, 0.88475, 
                                       0.885492, 0.886234, 0.886977, 0.887721, 0.888465, 0.88921, 0.889956, 
                                       0.890702, 0.891449, 0.892196, 0.892944, 0.893692, 0.894441, 0.895191, 
                                       0.895941, 0.896692, 0.897444, 0.898196, 0.898948, 0.899701, 0.900455, 
                                       0.901209, 0.901964, 0.90272, 0.903476, 0.904232, 0.904989, 0.905747, 
                                       0.906505, 0.907264, 0.908023, 0.908783, 0.909543, 0.910304, 0.911066, 
                                       0.911828, 0.91259, 0.913353, 0.914116, 0.914881, 0.915645, 0.91641, 
                                       0.917176, 0.917942, 0.918708, 0.919475, 0.920243, 0.921011, 0.921779, 
                                       0.922548, 0.923318, 0.924087, 0.924858, 0.925629, 0.9264, 0.927172, 
                                       0.927944, 0.928716, 0.929489, 0.930263, 0.931037, 0.931811, 0.932586, 
                                       0.933361, 0.934137, 0.934913, 0.935689, 0.936466, 0.937243, 0.938021, 
                                       0.938799, 0.939577, 0.940356, 0.941135, 0.941915, 0.942695, 0.943475, 
                                       0.944256, 0.945037, 0.945818, 0.946599, 0.947381, 0.948164, 0.948946, 
                                       0.949729, 0.950512, 0.951296, 0.95208, 0.952864, 0.953648, 0.954433, 
                                       0.955218, 0.956003, 0.956789, 0.957574, 0.95836, 0.959147, 0.959933, 
                                       0.96072, 0.961507, 0.962294, 0.963081, 0.963869, 0.964657, 0.965445, 
                                       0.966233, 0.967021, 0.96781, 0.968599, 0.969387, 0.970176, 0.970966, 
                                       0.971755, 0.972544, 0.973334, 0.974124, 0.974913, 0.975703, 0.976493, 
                                       0.977283, 0.978074, 0.978864, 0.979654, 0.980445, 0.981235, 0.982025, 
                                       0.982816, 0.983606, 0.984397, 0.985187, 0.985978, 0.986768, 0.987559, 
                                       0.988349, 0.98914, 0.98993, 0.990721, 0.991511, 0.992301, 0.993091, 
                                       0.993881, 0.994671, 0.995461, 0.996251, 0.99704, 0.99783, 0.998619)
      lines(RateofConsumption, TrappingProbabilityOriginal, type = "l", lty = MyLines[k], lwd = 2, col = "red")
    } else if ( k == 3){
      TrappingProbabilityCashTransfer <- c(0.69006, 0.69113, 0.692202, 0.693277, 0.694355, 0.695435, 0.696518, 
                                           0.697604, 0.698692, 0.699783, 0.700877, 0.701974, 0.703073, 0.704175, 
                                           0.70528, 0.706387, 0.707497, 0.70861, 0.709726, 0.710845, 0.711966, 
                                           0.71309, 0.714217, 0.715347, 0.716479, 0.717615, 0.718753, 0.719894, 
                                           0.721038, 0.722185, 0.723335, 0.724487, 0.725643, 0.726801, 0.727963, 
                                           0.729127, 0.730294, 0.731464, 0.732637, 0.733813, 0.734992, 0.736174, 
                                           0.737358, 0.738546, 0.739737, 0.740931, 0.742128, 0.743328, 0.744531, 
                                           0.745736, 0.746945, 0.748157, 0.749372, 0.750591, 0.751812, 0.753036, 
                                           0.754263, 0.755494, 0.756728, 0.757964, 0.759204, 0.760447, 0.761693, 
                                           0.762943, 0.764195, 0.765451, 0.76671, 0.767972, 0.769237, 0.770505, 
                                           0.771777, 0.773052, 0.77433, 0.775611, 0.776896, 0.778184, 0.779475, 
                                           0.78077, 0.782067, 0.783368, 0.784673, 0.78598, 0.787291, 0.788606, 
                                           0.789923, 0.791244, 0.792569, 0.793897, 0.795228, 0.796562, 0.7979, 
                                           0.799242, 0.800586, 0.801935, 0.803286, 0.804641, 0.806, 0.807362, 
                                           0.808728, 0.810097, 0.811469, 0.812845, 0.814225, 0.815608, 0.816994, 
                                           0.818384, 0.819778, 0.821175, 0.822576, 0.82398, 0.825388, 0.8268, 
                                           0.828215, 0.829633, 0.831056, 0.832482, 0.833911, 0.835345, 0.836782, 
                                           0.838222, 0.839667, 0.841115, 0.842566, 0.844022, 0.845481, 0.846944, 
                                           0.84841, 0.84988, 0.851354, 0.852832, 0.854314, 0.855799, 0.857288, 
                                           0.858781, 0.860278, 0.861779, 0.863283, 0.864791, 0.866303, 0.867819, 
                                           0.869339, 0.870863, 0.87239, 0.873922, 0.875457, 0.876996, 0.878539, 
                                           0.880086, 0.881637, 0.883192, 0.884751, 0.886314, 0.887881, 0.889451, 
                                           0.891026, 0.892605, 0.894188, 0.895774, 0.897365, 0.89896, 0.900559, 
                                           0.902162, 0.903769, 0.90538, 0.906995, 0.908614, 0.910237, 0.911865, 
                                           0.913496, 0.915132, 0.916771, 0.918415, 0.920063, 0.921715, 0.923371, 
                                           0.925032, 0.926696, 0.928365, 0.930038, 0.931715, 0.933396, 0.935082, 
                                           0.936771, 0.938465, 0.940163, 0.941866, 0.943572, 0.945283, 0.946998, 
                                           0.948718, 0.950441, 0.952169, 0.953901, 0.955638, 0.957379, 0.959124, 
                                           0.960873, 0.962627, 0.964385, 0.966147, 0.967914, 0.969685, 0.97146, 
                                           0.97324, 0.975024, 0.976813, 0.978605, 0.980403, 0.982204, 0.98401, 
                                           0.98582, 0.987635, 0.989454, 0.991278, 0.993106, 0.994938, 0.996775)
      lines(RateofConsumption, TrappingProbabilityCashTransfer, type = "l", lty = MyLines[k], lwd = 2, col = "blue")
      TrappingProbabilityOriginal <- c(0.699933, 0.700985, 0.702041, 0.703098, 0.704158, 0.705221, 0.706286, 
                                       0.707353, 0.708423, 0.709496, 0.710571, 0.711649, 0.712729, 0.713812, 
                                       0.714897, 0.715985, 0.717076, 0.718169, 0.719264, 0.720363, 0.721464, 
                                       0.722567, 0.723673, 0.724782, 0.725893, 0.727007, 0.728124, 0.729243, 
                                       0.730365, 0.73149, 0.732617, 0.733747, 0.73488, 0.736015, 0.737153, 
                                       0.738294, 0.739437, 0.740583, 0.741732, 0.742884, 0.744038, 0.745196, 
                                       0.746356, 0.747518, 0.748684, 0.749852, 0.751023, 0.752197, 0.753374, 
                                       0.754553, 0.755735, 0.75692, 0.758108, 0.759299, 0.760493, 0.761689, 
                                       0.762889, 0.764091, 0.765296, 0.766504, 0.767715, 0.768929, 0.770146, 
                                       0.771365, 0.772588, 0.773813, 0.775042, 0.776273, 0.777508, 0.778745, 
                                       0.779985, 0.781229, 0.782475, 0.783724, 0.784976, 0.786232, 0.78749, 
                                       0.788751, 0.790015, 0.791283, 0.792553, 0.793827, 0.795103, 0.796383, 
                                       0.797666, 0.798951, 0.80024, 0.801532, 0.802827, 0.804126, 0.805427, 
                                       0.806731, 0.808039, 0.80935, 0.810664, 0.811981, 0.813301, 0.814624, 
                                       0.815951, 0.817281, 0.818614, 0.81995, 0.821289, 0.822632, 0.823978, 
                                       0.825327, 0.826679, 0.828035, 0.829394, 0.830756, 0.832122, 0.83349, 
                                       0.834862, 0.836238, 0.837616, 0.838998, 0.840383, 0.841772, 0.843164, 
                                       0.844559, 0.845958, 0.84736, 0.848765, 0.850174, 0.851586, 0.853002, 
                                       0.854421, 0.855843, 0.857269, 0.858698, 0.860131, 0.861567, 0.863006, 
                                       0.864449, 0.865896, 0.867346, 0.868799, 0.870256, 0.871716, 0.87318, 
                                       0.874647, 0.876118, 0.877593, 0.87907, 0.880552, 0.882037, 0.883525, 
                                       0.885017, 0.886513, 0.888012, 0.889515, 0.891021, 0.892531, 0.894045, 
                                       0.895562, 0.897083, 0.898607, 0.900135, 0.901667, 0.903202, 0.904741, 
                                       0.906284, 0.90783, 0.90938, 0.910933, 0.912491, 0.914051, 0.915616, 
                                       0.917184, 0.918756, 0.920332, 0.921912, 0.923495, 0.925082, 0.926672, 
                                       0.928267, 0.929865, 0.931467, 0.933072, 0.934682, 0.936295, 0.937912, 
                                       0.939533, 0.941157, 0.942785, 0.944418, 0.946053, 0.947693, 0.949337, 
                                       0.950984, 0.952635, 0.95429, 0.955949, 0.957612, 0.959278, 0.960949, 
                                       0.962623, 0.964301, 0.965983, 0.967669, 0.969359, 0.971052, 0.97275, 
                                       0.974451, 0.976156, 0.977866, 0.979579, 0.981295, 0.983016, 0.984741, 
                                       0.98647, 0.988202, 0.989939, 0.991679, 0.993424, 0.995172, 0.996924)
      lines(RateofConsumption, TrappingProbabilityOriginal, type = "l", lty = MyLines[k], lwd = 2, col = "red")
    } else if ( k == 4){
      TrappingProbabilityCashTransfer <- c(0.639473, 0.640644, 0.641818, 0.642996, 0.644177, 0.645361, 0.646549, 
                                           0.64774, 0.648935, 0.650134, 0.651335, 0.652541, 0.653749, 0.654962, 
                                           0.656178, 0.657397, 0.65862, 0.659847, 0.661077, 0.662311, 0.663548, 
                                           0.664789, 0.666034, 0.667282, 0.668534, 0.66979, 0.67105, 0.672313, 
                                           0.67358, 0.67485, 0.676125, 0.677403, 0.678685, 0.679971, 0.68126, 
                                           0.682553, 0.683851, 0.685152, 0.686457, 0.687766, 0.689078, 0.690395, 
                                           0.691715, 0.69304, 0.694368, 0.695701, 0.697037, 0.698378, 0.699722, 
                                           0.701071, 0.702423, 0.70378, 0.70514, 0.706505, 0.707874, 0.709247, 
                                           0.710624, 0.712005, 0.71339, 0.71478, 0.716173, 0.717571, 0.718973, 
                                           0.72038, 0.72179, 0.723205, 0.724624, 0.726048, 0.727476, 0.728908, 
                                           0.730344, 0.731785, 0.73323, 0.73468, 0.736134, 0.737592, 0.739055, 
                                           0.740522, 0.741994, 0.74347, 0.744951, 0.746436, 0.747926, 0.749421, 
                                           0.750919, 0.752423, 0.753931, 0.755444, 0.756961, 0.758483, 0.76001, 
                                           0.761541, 0.763078, 0.764618, 0.766164, 0.767714, 0.769269, 0.770829, 
                                           0.772394, 0.773964, 0.775538, 0.777117, 0.778701, 0.78029, 0.781884, 
                                           0.783483, 0.785087, 0.786695, 0.788309, 0.789928, 0.791552, 0.79318, 
                                           0.794814, 0.796453, 0.798097, 0.799746, 0.8014, 0.80306, 0.804724, 
                                           0.806394, 0.808069, 0.809749, 0.811434, 0.813125, 0.814821, 0.816522, 
                                           0.818229, 0.81994, 0.821658, 0.82338, 0.825108, 0.826841, 0.82858, 
                                           0.830324, 0.832074, 0.833829, 0.83559, 0.837356, 0.839128, 0.840905, 
                                           0.842688, 0.844476, 0.84627, 0.84807, 0.849875, 0.851686, 0.853503, 
                                           0.855325, 0.857153, 0.858987, 0.860827, 0.862672, 0.864524, 0.866381, 
                                           0.868244, 0.870112, 0.871987, 0.873868, 0.875754, 0.877647, 0.879545, 
                                           0.88145, 0.88336, 0.885277, 0.887199, 0.889128, 0.891063, 0.893003, 
                                           0.89495, 0.896903, 0.898863, 0.900828, 0.9028, 0.904778, 0.906762, 
                                           0.908753, 0.910749, 0.912752, 0.914762, 0.916778, 0.9188, 0.920828, 
                                           0.922863, 0.924905, 0.926952, 0.929007, 0.931068, 0.933135, 0.935209, 
                                           0.937289, 0.939376, 0.94147, 0.94357, 0.945677, 0.947791, 0.949911, 
                                           0.952038, 0.954172, 0.956312, 0.95846, 0.960614, 0.962774, 0.964942, 
                                           0.967117, 0.969298, 0.971486, 0.973682, 0.975884, 0.978093, 0.980309, 
                                           0.982532, 0.984762, 0.986999, 0.989243, 0.991494, 0.993753, 0.996018)
      lines(RateofConsumption, TrappingProbabilityCashTransfer, type = "l", lty = MyLines[k], lwd = 2, col = "blue")
    }
    TrappingProbabilityOriginal <- c(0.648622, 0.64978, 0.650941, 0.652105, 0.653272, 0.654442, 0.655616, 
                                     0.656793, 0.657974, 0.659157, 0.660344, 0.661535, 0.662728, 0.663925, 
                                     0.665126, 0.66633, 0.667537, 0.668747, 0.669961, 0.671179, 0.6724, 
                                     0.673624, 0.674852, 0.676083, 0.677318, 0.678556, 0.679798, 0.681044, 
                                     0.682292, 0.683545, 0.684801, 0.68606, 0.687324, 0.68859, 0.689861, 
                                     0.691135, 0.692412, 0.693694, 0.694979, 0.696267, 0.69756, 0.698856, 
                                     0.700156, 0.701459, 0.702766, 0.704077, 0.705392, 0.706711, 0.708033, 
                                     0.709359, 0.710689, 0.712023, 0.71336, 0.714702, 0.716047, 0.717397, 
                                     0.71875, 0.720107, 0.721468, 0.722833, 0.724202, 0.725575, 0.726952, 
                                     0.728332, 0.729717, 0.731106, 0.732499, 0.733896, 0.735297, 0.736702, 
                                     0.738112, 0.739525, 0.740943, 0.742364, 0.74379, 0.74522, 0.746654, 
                                     0.748092, 0.749535, 0.750982, 0.752433, 0.753888, 0.755347, 0.756811, 
                                     0.758279, 0.759752, 0.761229, 0.76271, 0.764195, 0.765685, 0.767179, 
                                     0.768678, 0.770181, 0.771688, 0.7732, 0.774717, 0.776238, 0.777763, 
                                     0.779293, 0.780827, 0.782366, 0.78391, 0.785458, 0.787011, 0.788568, 
                                     0.79013, 0.791696, 0.793267, 0.794843, 0.796424, 0.798009, 0.799599, 
                                     0.801194, 0.802793, 0.804397, 0.806006, 0.80762, 0.809239, 0.810862, 
                                     0.812491, 0.814124, 0.815762, 0.817405, 0.819052, 0.820705, 0.822363, 
                                     0.824026, 0.825693, 0.827366, 0.829044, 0.830726, 0.832414, 0.834107, 
                                     0.835805, 0.837508, 0.839216, 0.840929, 0.842647, 0.844371, 0.846099, 
                                     0.847833, 0.849572, 0.851317, 0.853066, 0.854821, 0.856581, 0.858347, 
                                     0.860118, 0.861894, 0.863675, 0.865462, 0.867254, 0.869052, 0.870855, 
                                     0.872663, 0.874477, 0.876297, 0.878122, 0.879952, 0.881788, 0.88363, 
                                     0.885477, 0.88733, 0.889188, 0.891052, 0.892921, 0.894796, 0.896677, 
                                     0.898564, 0.900456, 0.902354, 0.904258, 0.906167, 0.908083, 0.910004, 
                                     0.911931, 0.913864, 0.915802, 0.917747, 0.919697, 0.921653, 0.923615, 
                                     0.925584, 0.927558, 0.929538, 0.931524, 0.933516, 0.935514, 0.937518, 
                                     0.939529, 0.941545, 0.943567, 0.945596, 0.947631, 0.949672, 0.951719, 
                                     0.953772, 0.955832, 0.957897, 0.959969, 0.962048, 0.964132, 0.966223, 
                                     0.96832, 0.970424, 0.972534, 0.97465, 0.976772, 0.978902, 0.981037, 
                                     0.983179, 0.985327, 0.987482, 0.989644, 0.991812, 0.993986, 0.996167)
    lines(RateofConsumption, TrappingProbabilityOriginal, type = "l", lty = MyLines[k], lwd = 2, col = "red")
  }
  dev.off()
}

SensitivityAnalysisforRateofConsumptionTrappingProbability(InitialCapitals = c(1.3, 1.7, 4, 6), file = "/Users/josemiguelflorescontro/Desktop")



