# setting work space
setwd('C://Users/Clancy/Desktop/brainGraph_tutorial/')
# import packages
library(data.table)
## get cores for multiple processes
suppressMessages(library(brainGraph))
OS <- .Platform$OS.type
if (OS == 'windows'){
  pacman::p_load(snow, doSNOW)
  num.cores <- as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS'))
  cl <- makeCluster(num.cores, type='SOCK')
  clusterExport(cl, 'sim.rand.graph.par')
  registerDoSNOW(cl)
} else {
  suppressMessages(library(doMC))
  registerDoMC(detectCores())
}
## original matrix inspect
## shape = [164, 164], range = [32.87, 0]
#sub001 = read.csv('sub-sub001_parcels_2009.csv', header = FALSE) 
#rownames(sub001) <- colnames(sub001)

###############################################################################
#    DTI Matrix Analysis
###############################################################################
options(bg.subject_id='ID', bg.group='group')
# step 0. data sorted
datadir = 'C://Users/Clancy/Desktop/Structure_Dataset_71/txtdata/'
"if (!file.exists(datadir)) dir.create(datadir)
CsvFiles = list.files(path='.', pattern = 'sub-sub', full.names = T)
left_cortices = seq(1,74,1); right_cortices = seq(90,163,1)
total_cortices = c(left_cortices, right_cortices)
for (i in CsvFiles){
  OldFileName = tail(strsplit(i, '/')[[1]], 1)
  InputData = read.csv(i, header = FALSE)[total_cortices, total_cortices]
  NewFileName = sub('.csv', '.txt', OldFileName)
  write.table(InputData, file = NewFileName, row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
}

# >* means match all. ^ for begin and $ for end
file.rename(list.files(path = '.', pattern = '^sub-sub.*.txt$'),
            paste0(datadir, list.files(path = '.', pattern = '^sub-sub.*.txt$')))"

#Step 1. import data and create connectivity matrices
grps = c('health', 'patient')
## 1.1 covariance
covars.original = fread(file.path('C://Users/Clancy/Desktop/Structure_Dataset_71/subjects_71_info.csv'))
covars.dti = covars.original[, -c('disease_month')]
covars.dti[, group := as.factor(group)]
covars.dti[, gender := as.factor(gender)]
setkey(covars.dti, ID) # set key for sorting by the attribution
matfiles <- 
  list(A = list.files(datadir, pattern='sub-sub', full.names = T),
       label = list.files(datadir, pattern='label.txt', full.names = T))
inds = lapply(grps, function(x) covars.dti[group == x, which = TRUE]) # index for groups
thresholds <- seq(0.01, 0.5, 0.01)
my.mats <- create_mats(matfiles$A, modality = 'dti', divisor = 'none',
                        mat.thresh = thresholds, inds = inds)
A.norm.sub <- my.mats$A.norm.sub
A.norm.mean <- my.mats$A.norm.mean

# Step 2. create graphs and calculate metrics
setkey(destrieux, index)
## setorder(destrieux, index)  is similar with setkey.
atlas = 'destrieux'
g <- g.group <- vector('list', length(thresholds))
for (j in seq_along(thresholds)) {
  g[[j]] <- make_brainGraphList(A.norm.sub[[j]], atlas, modality = 'dti', 
                                weighting = 'sld', threshold = thresholds[j],
                                weighted = TRUE, gnames = covars.dti$ID,
                                grpNames = covars.dti$group )
  g.group[[j]] <- make_brainGraphList(A.norm.mean[[j]], atlas, modality = 'dti',
                                  weighting = 'sld', threshold = thresholds[j],
                                  weighted = TRUE, grpNames = grps )
}

## Graph- and vertex-level measures
dt.V <- rbindlist(lapply(g, vertex_attr_dt))
dt.G <- rbindlist(lapply(g, graph_attr_dt))
dt.V.group <- rbindlist(lapply(g.group, vertex_attr_dt))
dt.G.group <- rbindlist(lapply(g.group, graph_attr_dt))
## G-
idvars <- c('modality', 'atlas', 'weighting', 'threshold', 'Group')
dt.G.long <- melt(dt.G, id.vars = c(idvars, 'Study.ID', 'density'))
dt.G.group.long <- melt(dt.G.group, id.vars = c(idvars, 'Study.ID', 'density'))

#### mean degree for the group-averaged graphs
dt.V.group[, .(mean.degree = mean(degree)), by = .(Group, threshold, lobe)]
#### t-test of nodal efficiency for Frontal lobe vertices
dt.V.group[lobe=='Frontal', .(p=t.test(E.nodal~Group)$p.value), by=threshold]
dt.V.group[lobe=='Parietal', .(p=t.test(E.nodal~Group)$p.value), by=threshold]
dt.V.group[lobe=='Temporal', .(p=t.test(E.nodal~Group)$p.value), by=threshold]
dt.V.group[lobe=='Occipital', .(p=t.test(E.nodal~Group)$p.value), by=threshold] #0.002
dt.V.group[lobe=='Insula', .(p=t.test(E.nodal~Group)$p.value), by=threshold]
dt.V.group[lobe=='Limbic', .(p=t.test(E.nodal~Group)$p.value), by=threshold]
#### plot
#exclude.vars <- c('assortativity.lobe.hemi', 'max.comp', 'diameter.wt', 
#                  'clique.num', 'num.tri') 
#plot_global(g, xvar='threshold', exclude=exclude.vars)

# Step 3. Group Analyses: GLM-based
### dummy coding example
X <- covars.dti[, 2:3]
X[, levels(group)] # set group to factor
model.matrix( ~group, data=X)[c(1:4, 21:24), ] # dummy coding
# first nodal efficiency for each subject
y <- sapply(g[[1]]$graphs, function(x) V(x)$E.nodal.wt[1]) 
Xy <- cbind(X,y)
summary(lm(y~group, data=Xy))$coefficients # GLM model
# beta_0 equals the mean of Group 1, and beta_1 equals the difference in Group 
# means (Group 2 - Group 1)
Xy[, mean(y), by=group] # mean value of each group
Xy[, mean(y), by=group][, diff(V1)] # difference value of two groups

### effect coding example
Xdes <- brainGraph_GLM_design(X, coding='effects', factorize=FALSE)#effect coding
Xdes[c(1:4, 21:24), ]
# contrasts
Cmat <- matrix(c(1, -1, 1, 1, 0, -2, 0, 2), nrow = 4, ncol = 2, byrow = TRUE)
rownames(Cmat) <- c(paste('Mean', grps),
                    'Control > Patient', 'Patient > Control')
# fit with GLM
fits <- fastLmBG(Xdes, as.matrix(y))
rbindlist(apply(fastLmBG_t(fits, Cmat), 3, as.data.table), idcol='Contrast')

### cell means example
model.matrix( ~group+0, data=X)[c(1:4, 21:24), ] # cell means
summary(lm(y~group+0, data=Xy))$coefficients

##############################################################################
## myself GLM
##############################################################################

g.glm <- g[[1]] # g[[1]] means first threshold, it can be other number.

## base information
BaseInformation <- covars.dti[, 2:4]
#con.mat <- matrix(c(0, 0, -2), nrow=1, dimnames=list('Control > Patient'))
CompareMatrix <- matrix(c(0, 0, -2, 0, 0, 2), nrow = 2, ncol = 3, byrow = TRUE)
rownames(CompareMatrix) <- c('Control > Patient', 'Patient > Control')
summary(brainGraph_GLM(g[[1]], measure = 'E.nodal.wt', covars = BaseInformation,
               coding='effects', mean.center=TRUE, contrasts = CompareMatrix, alt = 'two.sided', permute = TRUE))

## Two-group difference with continuous covariate interaction
BaseInformation <- covars.dti[, 2:4]
X <- brainGraph_GLM_design(BaseInformation, coding = 'cell.means', mean.center = TRUE, int=c('group', 'age'))
CompareMatrix <- matrix(c(0, 0, 1, -1), nrow=1, dimnames=list('Group X Age'))

summary(brainGraph_GLM(g[[1]], measure = 'E.nodal.wt', covars = BaseInformation,
                       X=X, contrasts = CompareMatrix, alt = 'two.sided'))

## Two-group difference with binary covariate interaction
BaseInformation <- covars.dti[, c(2,3,5)]
X <- brainGraph_GLM_design(BaseInformation, coding = 'cell.means', 
                           center.by = c('group', 'gender'), int = c('group', 'gender'))
CompareMatrix <- matrix(c(0, 0, 1, -1), nrow=1, dimnames=list('Group X Gender'))

summary(brainGraph_GLM(g[[1]], measure = 'E.nodal.wt', covars = BaseInformation,
                       X=X, contrasts = CompareMatrix, alt = 'two.sided'))

## GLM with permutation test
WK <- covars.dti[, 2:6]
X <- brainGraph_GLM_design(WK, coding = 'effects',factorize = TRUE, binarize = 'gender')
con.mat <- matrix(c(rep(0, 4), -2), nrow = 1, dimnames = list('Control > Patient'))
diffs.perm <- brainGraph_GLM(g[[1]], measure = 'E.nodal.wt', covars = WK,
                             contrasts = con.mat, X=X, alt='greater', 
                             permute = TRUE, part.method = 'guttman', long = TRUE)
summary(diffs.perm) # part.method=c('beckmann', 'guttman', 'ridgway') not everyone are suit your data

# Step 4. Visualize results
library(gridExtra)
grid.arrange(plot(diffs.perm, region='lG_Ins_lg_and_S_cent_ins', which=1:6)[[1]])

##############################################################################
## MTPC multiple threshold permutation correction
##############################################################################
# clear all above
covars <- covars.dti[, 2:6]
mtpcVars <- data.table(level=rep(c('graph', 'graph'), each=2), 
                       outcome=c('E.global.wt', 'Lp', 'Cp', 'E.local.wt'), 
                       alt='greater')
# Change H_A for 'Lp' 
mtpcVars[outcome == 'Lp', alt := 'less'] 
setkey(mtpcVars, level, outcome)
# Different number of permutations based on the level 
mtpcVars['graph', N := 1e3] 
# Generate permutation matrices using 'shuffleSet' from the 'permute' package 
mtpcPerms <- list( 
                  graph=shuffleSet(n=nrow(covars), nset=mtpcVars['graph', unique(N)]))
# Create the contrast matrix 
mtpcContrast <- matrix(c(0, 0, 0, 0, -2, 0, 0, 0, 1, 0),nrow=2,byrow=TRUE, 
                       dimnames=list(c('Control vs. Patient', 'BMI effect')))
mtpc.diffs.list <- sapply(mtpcVars[, unique(level)], function(x) NULL)

for (x in names(mtpc.diffs.list)) { # Loop across 'level' 
  # The 2nd-level is for each network metric (for the given level 'x') 
  mtpc.diffs.list[[x]] <- sapply(mtpcVars[x, unique(outcome)], function(x) NULL) 
  for (y in mtpcVars[x, outcome]) { 
    # Print some timing info in the terminal; optional 
    print(paste('Level:', x, '; Outcome:', y, ';', format(Sys.time(), '%H:%M:%S'))) 
    
    mtpc.diffs.list[[x]][[y]] <- 
      mtpc(g, thresholds, covars=covars, measure=y, contrasts=mtpcContrast, 
           con.type='t', level=x, N=mtpcVars[.(x, y), N], perms=mtpcPerms[[x]],
           part.method = 'guttman', alt=mtpcVars[.(x, y), alt]) 
  } 
}

# data result display
mtpc.diffs.sig.dt <- 
  rbindlist(lapply(mtpc.diffs.list, function(x) 
    rbindlist(lapply(x, function(y) 
      y$DT[A.mtpc > A.crit, .SD[1], by=region]))))


################################################################################
## NBS network-based analysis
###############################################################################
BaseInformation <- covars.dti[, 2:6]
X <- brainGraph_GLM_design(BaseInformation, coding = 'effects', binarize = 'gender')
con.mat <- matrix(c(rep(0, 4), 2), nrow = 1, dimnames = list('Control > Patient '))
res.nbs <- NBS(my.mats$A, covars = BaseInformation, contrasts = con.mat, X=X,
                  p.init = 0.001, N=1e3, alternative ='greater', long = TRUE,
                  part.method = 'guttman')
WkMatreix = res.nbs$T.mat
WkMatreix = unlist(WkMatreix)
summary(res.nbs)
atlas = 'destrieux'
g.nbs <- make_brainGraphList(res.nbs, atlas=atlas)
# visualization
plot(g.nbs[1],alpha=0.05,label=atlas,vertex.label='lobe') # display with lobe name

# inspect elements in result
igraph_result <- g.nbs$graphs$`Control > Patient `
V(igraph_result)$p.nbs # p value of each node
V(igraph_result)$name[ V(igraph_result)$p.nbs > 0] # specific nodes with p value
E(igraph_result) # edges
E(igraph_result)$weight # edge weights
E(igraph_result)$dist # edge distances

################################################################################
## Mediation analysis
#medVars[.('graph', x, y), on=c('level','outcome','mediator'), treat]
#medVars[.('graph', x, y), on=c('level','outcome','mediator'), outcome]
#medVars[.('graph', x, y), on=c('level','outcome','mediator'), N]
###############################################################################
medVars <- data.table(level=c(rep('vertex', 2)), 
                      outcome='PRI', 
                      mediator=c('E.nodal.wt', 'E.local.wt'), 
                      treat='group', 
                      interact=FALSE, 
                      N=rep(1e4, 2))
BaseInformation <- covars.dti[, c(2,3,4,5,6,8)] 
med.list.g <- sapply(medVars['vertex', on='level', unique(outcome)], function(x) NULL)

for (x in names(med.list.g)) { # x='sf-36'
  med.list.g[[x]] <- sapply(medVars[.('vertex', x), on=c('level', 'outcome'), unique(mediator)], function(x) NULL)
  for (y in names(med.list.g[[x]])) { # y='E.global.wt'
    print(paste('Outcome:', x, '; Mediator:', y, ';', format(Sys.time(), '%H:%M:%S')))
    med.list.g[[x]][[y]] <- vector('list', length = length(thresholds))
    for (z in seq_along(med.list.g[[x]][[y]])) { # z=1,2,,,50
      med.list.g[[x]][[y]][[z]] <- 
        brainGraph_mediate(g[[z]], BaseInformation, mediator=y,
                           level = 'vertex', # change it if you specify to graph
                           treat=medVars[.('vertex', x, y), on=c('level','outcome','mediator'), treat],
                           outcome=medVars[.('vertex', x, y), on=c('level','outcome','mediator'), outcome],
                           covar.names=c('age','gender','BMI'),
                           boot=T, boot.ci.type='perc', N=medVars[.('vertex', x, y), on=c('level','outcome','mediator'), N],
                           long=TRUE, binarize=c('gender'))
    }
  }
}
# you can only run single with 'graph' or 'vertex', or you can add valiable J to for-loop.
