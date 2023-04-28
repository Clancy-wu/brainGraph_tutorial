library(brainGraph)
library(data.table)
options(bg.subject_id='participant_id', bg.group='group')
data.l.file <- fread('aparc.a2009s_lh_thickness.csv')
data.r.file <- fread('aparc.a2009s_rh_thickness.csv')
data.l <- data.l.file[,2:75]; data.r <- data.r.file[, 2:75]
lhrh <- cbind(data.l, data.r)
# sort name
library(stringr)
oldname <- colnames(lhrh)
newname <- str_replace_all(oldname, c('h_'='', '&'='_and_', '_thickness'='', '-'='.'))
colnames(lhrh) <- newname
# begin
covars.all <- fread('subjects_71_info.csv')
covars.all[, gender := as.factor(gender)]
covars <- covars.all[,c('age', 'gender', 'group')]
covars.id <- paste0('sub-', covars.all$participant_id)
covars$participant_id <- covars.id
lhrh$participant_id <- covars.id

myResids <- get.resid(lhrh, covars, atlas = 'destrieux')
densities <- seq(0.01, 0.34, 0.01) 
corrs <- corr.matrix(myResids, densities=densities)

g <- lapply(seq_along(densities), function(x) 
  make_brainGraphList(corrs[x], modality='thickness'))
dt.G <- rbindlist(lapply(g, graph_attr_dt))  
dt.V <- rbindlist(lapply(g, vertex_attr_dt))

plot(myResids, region='lS_calcarine', cols = T)[[1]]

setkey(dt.V, density, group)
dt.V[density < 0.1, .SD[which.max(degree), .(region, degree)], 
     by=.(group, density)]





