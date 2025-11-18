
######################################
####### SHERC EGA WORKSHOP ###########
####### Heidi A. Wayment #############
####### November 18, 2025 ############
######################################
# Check your working directory
getwd()

#####################################
########  Packages Needed  ##########
#####################################

install.packages(c("psych","EGAnet","qgraph", "dplyr", "haven", "lavaan"))
#omit any of these if already installed

library(psych); data(bfi) #data are in the psych package
library(dplyr)
library(EGAnet)
library(qgraph)
library(haven)
library(lavaan)
library(ggplot2)

########################################
######## BFI Data Information ##########
########################################

#Name: bfi (available in psych package)
#Source: International Personality Item Pool (IPIP; Goldberg, 1999), collected 
#via the Synthetic Aperture Personality Assessment (SAPA) project.
#Contents: N = 2800
#25 self-report items measuring the Big Five personality traits 
#(Openness, Conscientiousness, Extraversion, Agreeableness, Neuroticism)
#Includes three demographic variables: sex, education, age

#############################
######## Load Data ##########
#############################

data(bfi)
# take a look
head(bfi)

#if you want to save the data as a csv file
write.csv(bfi,"bfi_dataset.csv", row.names = FALSE)

# create a dataset called bfi_25 of the 25 bfi items

bfi_25 <- bfi |> dplyr::select(A1:A5, C1:C5, E1:E5, N1:N5, O1:O5) |> na.omit()
      #removed missing data now N = 2436 with 25 items
head(bfi_25)
dim(bfi_25)

#############################
#######     UVA   ###########
#############################

#Use UVA function in EGAnet to check for redundant items based on weighted topological overlap (wTO).
#Items with redundancy (ω ≥ 0.25) are flagged, and the function removes the more redundant ones.
#This step ensures cleaner networks by avoiding items that are essentially duplicates.

#run the UVA function 
bfi_uva <- UVA(
  data = bfi_25,
  key = as.character(bfi.dictionary$Item[1:25])
)

bfi_uva #see result

bfi_uva$keep_remove  # shows which items were removed

head(bfi_uva$reduced_data, 3) #just show 3 lines

#create a new datafile to use

bfi_red <- bfi_uva$reduced_data #option to save reduced data 
head(bfi_red, 3) #the red is for reduced

#################################################
######## BASIC EGA on Reduced Dataset  ##########
#################################################

ega_gl_red <- EGAnet::EGA(
  data = bfi_red, 
  algorithm = "louvain", 
  model = "glasso",
  corr = "auto",
  consensus.iter = 1000,
  consensus.method = c("highest_modularity"),
  plot.EGA = TRUE,
  seed = 1234
)

ega_gl_red #note the TEFI is -24.989 using the glasso estimation method


ega_t_red <- EGAnet::EGA(
  data = bfi_red, 
  algorithm = "louvain", 
  model = "tmfg",
  corr = "auto",
  consensus.iter = 1000,
  consensus.method = c("highest_modularity"),
  plot.EGA = TRUE,
  seed = 1234
)

ega_t_red #note the TEFI is -22.53 using the tmfg estimation method

ega_gl_red$TEFI #as in the full dataset, we'll continue with glasso estimation
ega_t_red$TEFI

# 5 communities (dimensions) found in the Glasso Model
# Matches the Five-Factor Model of personality

#################################################
######## Check Stability of Solution   ##########
#################################################

#Use bootEGA() to assess robustness of dimensions with bootstrap replicates (default = 500).
#Reports:
# Frequency of dimension recovery
# Median number of dimensions with 95% CI
# Item and dimension stability measures
# Structural consistency (how often a dimension is replicated exactly)

bfi_boot_r <- bootEGA(
  data = bfi_uva$reduced_data,
  seed = 1
)
summary(bfi_boot_r)
dimensionStability(bfi_boot_r)

plot.ega_r_boot <- plot(bfi_boot_r, vsize = 10, 
                        label.size = 4, 
                        alpha = 0.5, 
                        edge.alpha = 0.5, bg = "transparent",
                        legend.names = c("Agreeableness", "Conscienciousness",
                                         "Extraversion", "Neuroticism", "Openness"),
                        color.palette = "polychrome")

plot.ega_r_boot



############################################
# Let's print a better figure of our selected model #
############################################

plot.ega_glasso <- plot(ega_gl_red, vsize = 10, 
                        label.size = 4, 
                        alpha = 0.5, 
                        edge.alpha = 0.5, bg = "transparent",
                        legend.names = c("Agreeableness", "Conscienciousness",
                                         "Extraversion", "Neuroticism", "Openness"),
                        color.palette = "polychrome") + 
  ggtitle("BFI Data N = 2436, Glasso Estimation Method,
  on reduced set of 24 items")

plot.ega_glasso

ggsave("Final_Glasso_Plot.png", plot = plot.ega_glasso, bg = "transparent")

####################################################
######   COMPUTE A CFA ON REDUCED DATA EGA    ######
####################################################

ega_cfa_red <- EGAnet::CFA(ega.obj = ega_gl_red, estimator = "WLSMV",
                           plot.CFA = TRUE,
                           data = bfi_red)

lavaan::fitMeasures(ega_cfa_red$fit, fit.measures = c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

semPlot::semPaths(ega_cfa_red$fit,
                  label.cex = 0.8, sizeLat = 8, sizeMan = 5, 
                  edge.label.cex = 0.6, minimum = 0.1,
                  sizeInt = 0.8, mar = c(1, 1, 1, 1), residuals = FALSE, 
                  intercepts = FALSE, thresholds = FALSE, layout = "spring",
                  "std", cut = 0.5)

title("EFA on Reduced Dataset EGA-glasso",
      adj = 1, line = 3, cex.main = 1,   font.main= 2, col.main= "darkgreen",
      cex.sub = 1, font.sub = 4, col.sub = "black")


####Following code are same analyses on the full N = 25 item bfi data

#########################################################################
######## BASIC EGA : Trying Glasso and TMFG Estimation Methods ##########
#########################################################################

ega_glasso <- EGAnet::EGA(
      data = bfi_25, 
      algorithm = "louvain", 
      model = "glasso",
      corr = "auto",
      consensus.iter = 1000,
      consensus.method = c("highest_modularity"),
      plot.EGA = TRUE,
      seed = 1234
)

ega_glasso #note the TEFI is -27.501 using the glasso estimation method


ega_tmfg <- EGAnet::EGA(
      data = bfi_25, 
      algorithm = "louvain", 
      model = "tmfg",
      corr = "auto",
      consensus.iter = 1000,
      consensus.method = c("highest_modularity"),
      plot.EGA = TRUE,
      seed = 1234
)

ega_tmfg #note the TEFI is -26.11 using the tmfg estimation method

ega_glasso$TEFI #we'll continue with glasso estimation
ega_tmfg$TEFI


############################################
# Let's print a better labelled EGA figure #
############################################

plot.ega_glasso <- plot(ega_glasso, vsize = 10, 
                          label.size = 4, 
                          alpha = 0.5, 
                          edge.alpha = 0.5, bg = "transparent",
                          legend.names = c("Agreeableness", "Conscienciousness",
                                           "Extraversion", "Neuroticism", "Openness"),
                          color.palette = "polychrome") + 
                          ggtitle("BFI Data N = 2436, Glasso Estimation Method")

plot.ega_glasso




Assuming 'p' is your ggplot object
ggsave("my_transparent_plot.png", plot = p, bg = "transparent")

################################################################################
###### HOW STABLE IS THE STRUCTURE and the ITEMS in the Structure? #############
################################################################################

BOOTega_glasso <- EGAnet::bootEGA(
      data = bfi_25, 
      cor = "cor_auto",
      uni.method = "louvain",
      iter = 600, # Number of replica samples to generate; can reduce number
      # resampling" for n random subsamples of the original data
      # parametric" for n synthetic samples from multivariate normal dist.
      type = "parametric", 
      # EGA Uses standard exploratory graph analysis
      # EGA.fit Uses total entropy fit index (tefi) to determine best fit of EGA
      # hierEGA Uses hierarchical exploratory graph analysis
      EGA.type = "EGA", 
      model = "glasso", 
      algorithm = "walktrap", # or "louvain" (better for unidimensional structures)
      # use "highest_modularity", "most_common", or "lowest_tefi"
      consensus.method = "highest_modularity", 
      typicalStructure = TRUE, # typical network of partial correlations
      plot.typicalStructure = TRUE, # returns a plot of the typical network
      ncores = 7,
      progress = TRUE ,
      summary.table	= TRUE,
      seed = 1234
)

summary(BOOTega_glasso)
plot(BOOTega_glasso)

#plot the bootstrapped figure

plot.bootEGA <- plot(BOOTega_glasso, vsize = 10, 
                        label.size = 4, 
                        alpha = 0.5, 
                        edge.alpha = 0.5,
                        legend.names = c("Agreeableness", "Conscienciousness",
                                         "Extraversion", "Neuroticism", "Openness"),
                        color.palette = "polychrome") + 
  ggtitle("BFI Data N = 2436, Bootstrap 600 samples Estimation")

plot.bootEGA


dim_stability <- EGAnet::dimensionStability(BOOTega_glasso)
dim_stability
dim_stability$dimension.stability
dim_stability$item.stability

# Create the plot and store it
bf1_stab_info <- dim_stability$item.stability$plot

# Add a title using labs() or ggtitle()
bf1_stab_info + labs(title = "Bootstrapped Item Stability Results")

####################################################
######   COMPUTE A CFA BASED ON THE RESULTS   ######
####################################################

ega_cfa_glasso <- EGAnet::CFA(ega.obj = ega_glasso, estimator = "WLSMV",
                           plot.CFA = TRUE,
                           data = items)

ega_cfa_glasso

# specifying estimator = "WLSMV" is the recommended approach for categorical or ordinal variables
# compute fit measures with lavaan package

lavaan::fitMeasures(ega_cfa_glasso$fit, fit.measures = c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

semPlot::semPaths(ega_cfa_glasso$fit,
                  label.cex = 0.8, sizeLat = 8, sizeMan = 5, 
                  edge.label.cex = 0.6, minimum = 0.1,
                  sizeInt = 0.8, mar = c(1, 1, 1, 1), residuals = FALSE, 
                  intercepts = FALSE, thresholds = FALSE, layout = "spring",
                  "std", cut = 0.5)

title("EFA based on EGA Glasso Results", sub = "Nov 18 2025",
      adj = 1, line = 3, cex.main = 1,   font.main= 2, col.main= "darkgreen",
      cex.sub = 1, font.sub = 4, col.sub = "black")


############################################
##Does a higher-order structure fit better?
###########################################

Hierega_gl <- EGAnet::hierEGA(
      data = bfi_red, 
      algorithm = "louvain", 
      model = "glasso",
      corr = "auto",
      consensus.iter = 1000,
      consensus.method = c("highest_modularity"),
      plot.EGA = TRUE,
      seed = 1234
)

Hierega_gl$plot.hierEGA
Hierega_gl$TEFI

##Comparing if lower or higher order a better fit
Hierega_gl$lower_order$TEFI
Hierega_gl$higher_order$TEFI

################################################################################
######    HOW STABLE IS THE HIGHER-ORDER STRUCTURE and the ITEMS   #############
################################################################################

BOOTega_H <- EGAnet::bootEGA(
      data = bfi_red, 
      cor = "cor_auto",
      uni.method = "louvain",
      iter = 600, # Number of replica samples to generate
      # resampling" for n random subsamples of the original data
      # parametric" for n synthetic samples from multivariate normal dist.
      type = "parametric", 
      # EGA Uses standard exploratory graph analysis
      # EGA.fit Uses total entropy fit index (tefi) to determine best fit of EGA
      # hierEGA Uses hierarchical exploratory graph analysis
      EGA.type = "hierEGA", 
      model = "glasso", 
      algorithm = "walktrap", # or "louvain" (better for unidimensional structures)
      # use "highest_modularity", "most_common", or "lowest_tefi"
      consensus.method = "highest_modularity", 
      typicalStructure = TRUE, # typical network of partial correlations
      plot.typicalStructure = TRUE, # returns a plot of the typical network
      ncores = 7,
      progress = TRUE ,
      summary.table	= TRUE,
      seed = 1234
)

summary(BOOTega_H)


dim_stability_BH <- EGAnet::dimensionStability(BOOTega_H) #the figure
dim_stability_BH

BOOTega_H$EGA #suggests the lower order has a better fit.

#################
### Citations ###
#################

## Portions of this EGAnet workflow were adapted from the EGAnet documentation 
## developed by Hudson Golino and Alexander P. Christensen. 
## The EGAnet R package provides open resources for 
## applying Exploratory Graph Analysis (EGA) and related techniques.
## Golino, H. & Epskamp, S. (2017). Exploratory Graph Analysis: A new approach for
## estimating the number of dimensions in psychological research. PLOS ONE, 12(6), 
## e0174035. https://doi.org/10.1371/journal.pone.0174035

## Golino, H., Moulder, R., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D.,
## Thiyagarajan, J. A., & Martinez-Molina, A. (2020). Entropy fit index: 
##  A new fit measure for assessing the structure and dimensionality of multiple 
## latent variables. Multivariate Behavioral Research, 55(5), 556–583. 
## https://doi.org/10.1080/00273171.2020.1779642

## EGAnet documentation: https://r-ega.net
