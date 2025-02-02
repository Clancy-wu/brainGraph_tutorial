# 0_load_packages ---------------------------------------------------------
library(brainGraph) # compute network properties
# remotes::install_version("igraph", version = "1.6.0")
library(igraph) # 1.6.0
packageVersion("igraph") # 1.6.0
library(data.table) # data.table
library(parallel)
library(doMC)
registerDoMC(30)

options(bg.subject_id='participant_id', bg.group='all_group')
grps = c('health_before', 'patient_before', 'health_after', 'patient_after')
# 1_network_construction --------------------------------------------------

## Network contains positive & negative & absolute
densities <- seq(0.25, 0.30, 0.05)

## create new atlas
my_data = fread('brainnetome_surface.csv')
my_data$x.mni <- as.numeric(my_data$x.mni)
my_data$y.mni <- as.numeric(my_data$y.mni)
my_data$z.mni <- as.numeric(my_data$z.mni)
my_data$lobe <- as.factor(my_data$lobe)
my_data$hemi <- as.factor(my_data$hemi)
my_data$gyrus <- as.factor(my_data$gyrus)
my_data$Yeo_7network <- as.factor(my_data$Yeo_7network)
my_data$Yeo_17network <- as.factor(my_data$Yeo_17network)
setkey(my_data, index)
setindex(my_data, name)
my_data
brainnetome_surface = my_data
## save atlas
save(brainnetome_surface, file='brainnetome_surface.rda')
## load atlas
load('brainnetome_surface.rda')

## run with new atlas
for (i in seq_along(densities)){
  gw.sub[[i]] <- make_brainGraphList(my.mats$A.norm.sub[[i]], 'brainnetome_surface' , level='subject',
                                     modality = 'fmri',threshold = densities[i],
                                     weighted = TRUE, gnames = covars.all$participant_id,
                                     grpNames = covars.all$all_group )
}

## end. author@kangwu
