## --Diversity of visual inputs to Kenyon cells in the Drosophila mushroom body
## --Ganguly, Heckman, et al., 2024 

# The lines of code below use packages and functions originally developed by Bates et al., eLife 2020. 

#We show how we implemented them in order to access and visualize data in our figures

library(fafbseg)
library(dplyr)   

#--Inputs to KCs

#The 630 version data release (June 2023) was used.
download_flywire_release_data(version=630L)
flywire_connectome_data_version(set=630)

#Obtain data frames of all inputs to KCg-d and KCab-p. 
#Setting threshold to 4 returns all neurons connected by 5 or more synapses
#Use super_class designations in the resulting data frame to search for 'visual_projection' inputs
KCg_d_Inputs = flywire_partner_summary2("KCg-d", partners = 'in', threshold = 4)
KCg_d_Inputs

KCab_p_Inputs = flywire_partner_summary2("KCab-p", partners = 'in', threshold = 4)
KCab_p_Inputs



#--Adjacency Matrices for Fig 5 Supp 1A-C

#Neuroglancer scenes are used as inputs for the adjacency matrix function. 

#Left KCg-m's that get uniglomerular PN input = https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/5478112169033728 
#olfPNs with Input to Left KCgms = https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/5898988563726336
#OSNs with Input to the olfPNs above = https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/5545626437681152 
#Left KCg-d's (minus the ones that don't get visual input) = https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/4844894168285184
#VPNs and LVINs with direct input to left KCg-d's = https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/6433277864837120

VPNsLVINs_KCgd = flywire_adjacency_matrix(
  rootids = NULL,
  inputids = "https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/6433277864837120",
  outputids = "https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/4844894168285184",
  sparse = FALSE,
  Verbose = T
)

#Binarizes the matrix values so that connections with 5 or more synapses are set to 1. This does impact the hierarchical clustering in heatmap()
VPNsLVINs_KCgd_binarized = as.matrix((VPNsLVINs_KCgd>4)+0)

#Visualize binarized connections in a heatmap. All connections are set to black. 
heatmap(VPNsLVINs_KCgd_binarized, scale='none', keep.dendro = TRUE, col=c("white", "black"), labRow = FALSE, labCol = FALSE)


