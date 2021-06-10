library(Spectre)
library(parallel)

source("/project/chronoclust/tracksom_eval/eval_GP/FlowSOM-tracking/experiment_materials/evaluation/calculate_ari.R")
source("/project/chronoclust/tracksom_eval/eval_GP/FlowSOM-tracking/experiment_materials/evaluation/calculate_tracking_accuracy.R")
source("/project/chronoclust/tracksom_eval/eval_GP/FlowSOM-tracking/experiment_materials/evaluation/map_clusters.R")

# obtain data tables for each solution
gating.dir <- "/project/chronoclust/tracksom_eval/evaluation/Samples-with-noise/"
tracksom.dir <- "/project/chronoclust/tracksom_eval/tracksom_with_merge/clustering/"

args <- commandArgs(trailingOnly=TRUE)

option <- as.numeric(args[1])
tracksom.dir.options <- paste0(tracksom.dir, 'option_',option,'/')
output.dir <- "/project/chronoclust/tracksom_eval/eval_GP/with_merge/20210528/"
dir.create(output.dir, recursive=TRUE)

gating.list <- Spectre::read.files(file.loc = gating.dir,
                                   file.type = '.csv',
                                   do.embed.file.names = TRUE)
transitions <- read.csv("/project/chronoclust/tracksom_eval/evaluation/transition_rules.csv")

# repeat for tracksom options

top.scores <- list()
all.scores <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), c('parameter', 'ari_1', 'ari_2',
                                                                  'ari_3', 'ari_4', 'ari_5', 'ari_6', 'ari_7', 'ari_8', 'ari_mean'))
for (i in 0:99) {
  print(i)
  solution.dir <- paste0(tracksom.dir.options, i)
  solution.list <- Spectre::read.files(file.loc = solution.dir,
                                       file.type = '.csv',
                                       do.embed.file.names = TRUE)
  for (j in 1:8) {
    solution.list[[j]] <- data.frame(solution.list[[j]], gating.list[[j]]$PopName)
  }
  solution.list <- Spectre::do.merge.files(solution.list)
  names(solution.list)[length(names(solution.list))] <- 'PopName'
  solution.list <- solution.list[!(solution.list$PopName == 'None')]
  ari <- suppressMessages(data.frame(calculate_ari(dat = solution.list,
                                                   timepoints = c(1:8), 
                                                   timepoint_col = 'FileNo', 
                                                   cluster_col = 'TrackSOM_metacluster_lineage_tracking', 
                                                   population_col = 'PopName')))
  
  colnames(ari) <- c('ari_1', 'ari_2', 'ari_3', 'ari_4', 'ari_5', 'ari_6', 'ari_7', 'ari_8', 'ari_mean')
  acc <- suppressMessages(calculate_tracking_accuracy(dat = solution.list, 
                                                      timepoints = c(1:8), 
                                                      timepoint_col = 'FileNo', 
                                                      cluster_lineage_col = 'TrackSOM_metacluster_lineage_tracking',
                                                      population_col = 'PopName',
                                                      rule=transitions))
  mapping <- suppressMessages(map_clusters(dat = solution.list,
                                           timepoint_col = 'FileNo',
                                           timepoints = c(1:8),
                                           cluster_col = 'TrackSOM_metacluster_lineage_tracking', 
                                           population_col = 'PopName'))
  mapping_out_dir <- paste0(output.dir, "/mappings_", option)
  dir.create(mapping_out_dir, recursive=TRUE)
  setwd(mapping_out_dir)
  Spectre::write.files(mapping, i)
  
  # write out the legal and non legal transitions
  legal_transition_out_dir <- paste0(output.dir, "/legal_trans_", option)
  dir.create(legal_transition_out_dir, recursive=TRUE)
  setwd(legal_transition_out_dir)
  Spectre::write.files(acc$legal_transitions, i)
  
  # write out the legal and non legal transitions
  illegal_transition_out_dir <- paste0(output.dir, "/illegal_trans_", option)
  dir.create(illegal_transition_out_dir, recursive=TRUE)
  setwd(illegal_transition_out_dir)
  Spectre::write.files(acc$illegal_transitions, i)
  
  all.scores <- rbind(all.scores, c(parameter = i, ari, tracking = acc$tracking_accuracy))
}

setwd(output.dir)
Spectre::write.files(all.scores, paste0('all_scores_tracksom_option', option))
