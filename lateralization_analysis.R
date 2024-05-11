# 0_load_packages ---------------------------------------------------------
library(brainGraph) # compute network properties
# remotes::install_version("igraph", version = "1.6.0")
library(igraph) # 1.6.0
packageVersion("igraph") # 1.6.0
library(data.table) # data.table
library(pracma)     # compute AUC
library(freesurferformats)
library(parallel)
library(stringr)
library(doMC)
registerDoMC(detectCores())

options(bg.subject_id='participant_id', bg.group='group')
grps = c('health', 'patient')
# 1_network_construction --------------------------------------------------

## This is processed in python script. For functional network, the pearson 
## correlation between two nodes is defined as the edge. For structural network,
## the sfit2_invnodevol_count-sum between two nodes is defined as the edge.
## The diagonal in functional network is 0 and in structural network is 0.

## CFS
## Functional network: NetFunc / CFS / PearCor_sub-subxxx.txt
## Fiber network: NetFiber / CFS / Connectome_sub-subxxx.txt

## large population
## Functional network: NetFunc / ISYB / PearCor_sub-xxx.txt
## Fiber network: NetFiber / ISYB / Connectome_sub-xxx.txt

## split network to positive & negative & absolute
## NetResults/CFS /subxxx 
##                      /abs_fiber_xxx, abs_func_xxx, neg_fiber_xxx, neg_func_xxx, pos_fiber_xxx, pos_func_xxx
##           /ISYB /subxxx

## Script is: Atlas / atlas_transform.py

# 2_network_properties_compute --------------------------------------------
## network similarity

## Network contains positive & negative & absolute
densities <- seq(0.01, 0.34, 0.01)
atlas = 'brainnetome'
## --------------------------------------------------------------------------
save2dataframe <- function(data, file_path_suffix){
  vertex_df <- rbindlist(lapply(data, vertex_attr_dt))
  graph_df <- rbindlist(lapply(data, graph_attr_dt))
  fwrite(vertex_df, paste0(file_path_suffix, '_vertex.csv'), sep=',', col.names=TRUE)
  fwrite(graph_df, paste0(file_path_suffix, '_graph.csv'), sep=',', col.names=TRUE)
}
#### CFS: functional network graph weighted subject level
covars.all <- fread('ClinicalMeasurements/CFS_subinfo.csv') # subs=79
inds = lapply(grps, function(x) covars.all[group == x, which = TRUE]) # 1=health=39, 2=patient=40
matfiles <- paste0('NetResults/CFS/', covars.all$participant_id, '/abs_func_', covars.all$participant_id, '.txt') 
my.mats <- create_mats(matfiles, modality = 'fmri',threshold.by = 'density',
                       mat.thresh = densities, inds = inds)
gw.sub <- vector('list', length(densities)) # ws: weighted subject
for (i in seq_along(densities)){
  gw.sub[[i]] <- make_brainGraphList(my.mats$A.norm.sub[[i]], atlas, level='subject',
                                     modality = 'fmri',threshold = densities[i],
                                     weighted = TRUE, gnames = covars.all$participant_id,
                                     grpNames = covars.all$group )
  }
saveRDS(gw.sub, file=file.path('BrainGraphRDS/', 'CFS_func_weighted_subject_abs.rds'), compress = 'xz')
save2dataframe(gw.sub, file=file.path('BrainGraphRDS/', 'CFS_func_weighted_subject_abs'))

#### CFS: fiber network graph weighted subject level
matfiles <- paste0('NetResults/CFS/', covars.all$participant_id, '/abs_fiber_', covars.all$participant_id, '.txt') 
my.mats <- create_mats(matfiles, modality = 'fmri',threshold.by = 'density',
                       mat.thresh = densities, inds = inds)
gw.sub <- vector('list', length(densities)) # ws: weighted subject
for (i in seq_along(densities)){
  gw.sub[[i]] <- make_brainGraphList(my.mats$A.norm.sub[[i]], atlas, level='subject',
                                     modality = 'fmri',threshold = densities[i],
                                     weighted = TRUE, gnames = covars.all$participant_id,
                                     grpNames = covars.all$group )
}
saveRDS(gw.sub, file=file.path('BrainGraphRDS/', 'CFS_fiber_weighted_subject_abs.rds'), compress = 'xz')
save2dataframe(gw.sub, file=file.path('BrainGraphRDS/', 'CFS_fiber_weighted_subject_abs'))

## --------------------------------------------------------------------------
#### ISYB: functional network graph weighted subject level
covars.all <- fread('ClinicalMeasurements/ISYB_subinfo.csv') # subs=215
inds = lapply(grps, function(x) covars.all[group == x, which = TRUE]) # 1=health=215, 2=patient=0
matfiles <- paste0('NetResults/ISYB/', covars.all$participant_id, '/abs_func_', covars.all$participant_id, '.txt') 
my.mats <- create_mats(matfiles, modality = 'fmri',threshold.by = 'density',
                       mat.thresh = densities, inds = inds)
gw.sub <- vector('list', length(densities)) # ws: weighted subject
for (i in seq_along(densities)){
  gw.sub[[i]] <- make_brainGraphList(my.mats$A.norm.sub[[i]], atlas, level='subject',
                                     modality = 'fmri',threshold = densities[i],
                                     weighted = TRUE, gnames = covars.all$participant_id,
                                     grpNames = covars.all$group )
}
saveRDS(gw.sub, file=file.path('BrainGraphRDS/', 'ISYB_func_weighted_subject_abs.rds'), compress = 'xz')
save2dataframe(gw.sub, file=file.path('BrainGraphRDS/', 'ISYB_func_weighted_subject_abs'))

#### ISYB: fiber network graph weighted subject level
matfiles <- paste0('NetResults/ISYB/', covars.all$participant_id, '/abs_fiber_', covars.all$participant_id, '.txt') 
my.mats <- create_mats(matfiles, modality = 'fmri',threshold.by = 'density',
                       mat.thresh = densities, inds = inds)
gw.sub <- vector('list', length(densities)) # ws: weighted subject
for (i in seq_along(densities)){
  gw.sub[[i]] <- make_brainGraphList(my.mats$A.norm.sub[[i]], atlas, level='subject',
                                     modality = 'fmri',threshold = densities[i],
                                     weighted = TRUE, gnames = covars.all$participant_id,
                                     grpNames = covars.all$group )
}
saveRDS(gw.sub, file=file.path('BrainGraphRDS/', 'ISYB_fiber_weighted_subject_abs.rds'), compress = 'xz')
save2dataframe(gw.sub, file=file.path('BrainGraphRDS/', 'ISYB_fiber_weighted_subject_abs'))

# 3_structural_lateralization ---------------------------------------------





# 4_functional_compensation -----------------------------------------------




