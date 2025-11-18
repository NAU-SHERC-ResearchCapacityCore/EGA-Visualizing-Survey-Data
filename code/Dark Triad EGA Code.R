
# Packages
library(EGAnet)    # EGA and bootEGA
library(dplyr)
library(stringr)

####################################
####     Dark Triad Example     ####
####################################

#Read in your data (in following example the imported data called drk_triad)

#keeping items only
drk_triad <- drk_triad |> dplyr::select(M1:P9) |> na.omit()
#removed missing data now N = 2436 with 25 items
head(drk_triad)
dim(drk_triad)

# Reverse-code the following items
drk_triad$N2_R <-  6 - drk_triad$N2
drk_triad$N6_R <-  6 - drk_triad$N6
drk_triad$N8_R <-  6 - drk_triad$N8
drk_triad$P2_R <-  6 - drk_triad$P2
drk_triad$P7_R <-  6 - drk_triad$P7

head(drk_triad)

drk_triad_n <- drk_triad |> dplyr::select(M1:M9, N1, N3:N5, N7, N9:P1, P3:P6,
                                        P8:P7_R) |> na.omit()

head(drk_triad_n)

#########################################################################
######## BASIC EGA : Trying Glasso and TMFG Estimation Methods ##########
#########################################################################

ega_glasso_dt <- EGAnet::EGA(
  data = drk_triad_n, 
  algorithm = "louvain", 
  model = "glasso",
  corr = "auto",
  consensus.iter = 1000,
  consensus.method = c("highest_modularity"),
  plot.EGA = TRUE,
  seed = 1234
)

ega_glasso_dt #note the TEFI is -24.063 using the glasso estimation method 4 factors


ega_tmfg_dt <- EGAnet::EGA(
  data = drk_triad_n, 
  algorithm = "louvain", 
  model = "tmfg",
  corr = "auto",
  consensus.iter = 1000,
  consensus.method = c("highest_modularity"),
  plot.EGA = TRUE,
  seed = 1234
)

ega_tmfg_dt #note the TEFI is -19.739 using the tmfg estimation method 4 factors

ega_glasso_dt$TEFI #we'll continue with glasso estimation
ega_tmfg_dt$TEFI


#############################
#######     UVA   ###########
#############################


dark_uva <- UVA(
  data = drk_triad_n)
dark_uva$keep_remove  # shows which items were removed

head(dark_uva$reduced_data, 3) #just show 3 lines

dark_uva #only M1 and M7 have wT0 = .27, keep M1 (its not wise to tell your secrets)
# remove M7 (there are things you should hide)

##################


BOOTega_glasso_dt <- EGAnet::bootEGA(
  data = drk_triad_n, 
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

summary(BOOTega_glasso_dt)
plot(BOOTega_glasso_dt)

#plot the bootstrapped figure

library(ggplot2)

plot.bootEGA_dt <- plot(BOOTega_glasso_dt, vsize = 10, 
                        label.size = 4, 
                        alpha = 0.5, 
                        edge.alpha = 0.5,
                        legend.names = c("Factor", "Factor 2",
                                         "Factor 3"),
                        color.palette = "polychrome") + 
  ggtitle("Dark Triad Data N = 18,193, Bootstrap 600 samples Estimation")

plot.bootEGA_dt


dim_stability <- EGAnet::dimensionStability(BOOTega_glasso_dt)
dim_stability
dim_stability$dimension.stability
dim_stability$item.stability

#identifies some less stable items M3 and M4, N5 and N9

