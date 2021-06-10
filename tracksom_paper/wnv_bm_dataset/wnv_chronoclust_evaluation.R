# CHANGE: gating.dir, chronoclust.dir, transitions, mapping directory
# COPY/MOVE: best performing parameter sets, csv results of all sets after execution

library(Spectre)
library(dplyr)
library(parallel)

source("/project/chronoclust/tracksom_eval/eval_GP/FlowSOM-tracking/evaluation/calculate_ari.R")
source("/project/chronoclust/tracksom_eval/eval_GP/FlowSOM-tracking/evaluation/calculate_tracking_accuracy.R")
source("/project/chronoclust/tracksom_eval/eval_GP/FlowSOM-tracking/evaluation/map_clusters.R")

 # sample 100 solutions for chronoclust 
set.seed(42) 
numbers <- c(1:59, 61:200) 
x <- sample(numbers, 100)

gating.dir <- "/project/chronoclust/tracksom_eval/evaluation/Samples-with-noise/"
chronoclust.dir <- "/project/chronoclust/givanna/tracksom/chronoclust_eval/clustering/"
output.dir <- "/project/chronoclust/tracksom_eval/eval_GP/eval_20201021/chronoclust"
dir.create(output.dir, recursive=TRUE)

gating.list <- Spectre::read.files(file.loc = gating.dir,
                                   file.type = '.csv',
                                   do.embed.file.names = TRUE)

transitions <- read.csv("/project/chronoclust/tracksom_eval/evaluation/transition_rules.csv")

top.scores <- list()
all.scores <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), c('parameter', 'ari_1', 'ari_2',
              'ari_3', 'ari_4', 'ari_5', 'ari_6', 'ari_7', 'ari_8', 'ari_mean', 'tracking'))
for (i in x) {
  print(i)
  solution.dir <- paste0(chronoclust.dir, i)
  solution.list <- Spectre::read.files(file.loc = solution.dir,
                                     file.type = '.csv',
                                     do.embed.file.names = TRUE)
  solution.list <- solution.list[1:8]
  result.dat <- fread(paste0(solution.dir, '/result.csv'))
  for (j in 1:8) {
    solution.list[[j]] <- data.frame(solution.list[[j]], gating.list[[j]]$PopName)
    result.sub <- result.dat[result.dat$timepoint == j - 1,]
    solution.list[[j]] <- left_join(x = solution.list[[j]], y = result.sub[, c("tracking_by_lineage", "tracking_by_association")], 
                      by = c("cluster_id" = "tracking_by_lineage"))
  }
  solution.list <- Spectre::do.merge.files(solution.list)
  solution.list$tracking_by_association[is.na(solution.list$tracking_by_association)] <- "None"
  names(solution.list)[length(names(solution.list)) - 1] <- 'PopName'
  names(solution.list)[length(names(solution.list))] <- 'tracking_by_association'
  solution.list <- solution.list[!(solution.list$PopName == 'None')]
  ari <- suppressMessages(data.frame(calculate_ari(dat = solution.list,
                                timepoints = c(1:8), 
                                timepoint_col = 'FileNo', 
                                cluster_col = 'cluster_id', 
                                population_col = 'PopName')))
  
  colnames(ari) <- c('ari_1', 'ari_2', 'ari_3', 'ari_4', 'ari_5', 'ari_6', 'ari_7', 'ari_8', 'ari_mean')
  acc <- suppressMessages(calculate_tracking_accuracy(dat = solution.list, 
                                     timepoints = c(1:8), 
                                     timepoint_col = 'FileNo', 
                                     cluster_lineage_col = 'cluster_id',
                                     population_col = 'PopName',
                                     rule=transitions))
  mapping <- suppressMessages(map_clusters(dat = solution.list,
                                           timepoints = c(1:8),
                          timepoint_col = 'FileNo',
                          cluster_col = 'cluster_id', 
                          population_col = 'PopName'))
  mapping_out_dir <- paste0(output.dir, "/mappings")
  dir.create(mapping_out_dir, recursive=TRUE)
  setwd(mapping_out_dir)
  
  Spectre::write.files(mapping, i)
  
  # write out the legal and non legal transitions
  legal_transition_out_dir <- paste0(output.dir, "/legal_trans")
  dir.create(legal_transition_out_dir, recursive=TRUE)
  setwd(legal_transition_out_dir)
  Spectre::write.files(acc$legal_transitions, i)
  
  # write out the legal and non legal transitions
  illegal_transition_out_dir <- paste0(output.dir, "/illegal_trans")
  dir.create(illegal_transition_out_dir, recursive=TRUE)
  setwd(illegal_transition_out_dir)
  Spectre::write.files(acc$illegal_transitions, i)
  
  all.scores <- rbind(all.scores, c(parameter = i, ari, tracking = acc$tracking_accuracy))
}
setwd(output.dir)
Spectre::write.files(all.scores, 'all_scores_chronoclust')
