library(Spectre)

source('~/Documents/GitHub/FlowSOM-tracking/experiment_materials/evaluation/calculate_ari.R')
source('~/Documents/GitHub/FlowSOM-tracking/experiment_materials/evaluation/map_clusters.R')
source('~/Documents/GitHub/FlowSOM-tracking/experiment_materials/evaluation/calculate_tracking_accuracy.R')

# obtain data tables for each solution
gating.dir <- "~/Documents/GitHub/FlowSOM-tracking/experiment_materials/synthetic_dataset/dataset/csv_files_withLabel/"
tracksom.dir <- "~/Dropbox (Sydney Uni)/tracksom/synthetic_dataset/2021/synthetic_run_20210126/"
tracksom.dir.options <- c(paste0(tracksom.dir, 'option_1/'), paste0(tracksom.dir, 'option_2/'), paste0(tracksom.dir, 'option_3/'))
output.dir <- "~/Documents/GitHub/FlowSOM-tracking/experiment_materials/synthetic_dataset/result_qualities/"

gating.list <- Spectre::read.files(file.loc = gating.dir,
                                   file.type = '.csv',
                                   do.embed.file.names = TRUE)


transitions <- read.csv('~/Documents/GitHub/FlowSOM-tracking/experiment_materials/synthetic_dataset/dataset/transition_rules.csv')
# repeat for tracksom options
# 1 for AA, 2 for PI, and 3 for PV
option <- 1
setwd(tracksom.dir)
for (dir in tracksom.dir.options) {
  top.scores <- list()
  all.scores <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c('parameter', 'ari_1', 'ari_2',
								'ari_3', 'ari_4', 'ari_5', 'ari_mean'))
  for (i in 0:99) {
    message(paste("Evaluating option", option, "parameter", i))
    solution.dir <- paste0(dir, i)
    valid = TRUE
    if (file.exists(paste0(solution.dir, "/Result_FileName_synthetic_d0.csv"))) {
      solution.list <- Spectre::read.files(file.loc = solution.dir,
                                       file.type = '.csv',
                                       do.embed.file.names = TRUE)
      for (j in 1:5) {
      	if (nrow(solution.list[[j]]) == nrow(gating.list[[j]])) {
      	  solution.list[[j]] <- data.frame(solution.list[[j]], gating.list[[j]]$PopName)
      	} else {
      	  print(i)
      	  valid = FALSE
      	  break
	      }
      }
      if (valid) {
        solution.list <- Spectre::do.merge.files(solution.list)
        names(solution.list)[length(names(solution.list))] <- 'PopName'
        #solution.list <- solution.list[!(solution.list$PopName == 'Noise')]
        ari <- data.frame(calculate_ari(dat = solution.list,
    				  timepoints = c(1:5), 
    				  timepoint_col = 'FileNo', 
    				  cluster_col = 'TrackSOM_metacluster_lineage_tracking',
    				  population_col = 'PopName'))
      
        colnames(ari) <- c('ari_1', 'ari_2', 'ari_3', 'ari_4', 'ari_5', 'ari_mean')
        acc <- calculate_tracking_accuracy(dat = solution.list, 
                                           timepoints = c(1:5), 
                                           timepoint_col = 'FileNo', 
                                           cluster_lineage_col = 'TrackSOM_metacluster_lineage_tracking',
                                           population_col = 'PopName',
                                           rule=transitions)
        # write out mapping
        mapping <- map_clusters(dat = solution.list,
                              timepoint_col = 'FileNo',
                              timepoints = c(1:5),
                              cluster_col = 'TrackSOM_metacluster_lineage_tracking', 
                              population_col = 'PopName')
        mapping_out_dir <- paste0(output.dir, "/fine/mappings_", option)
        dir.create(mapping_out_dir, recursive=TRUE)
        setwd(mapping_out_dir)
        Spectre::write.files(mapping, i)
        
        # write out the legal and non legal transitions
        legal_transition_out_dir <- paste0(output.dir, "/fine/legal_trans_", option)
        dir.create(legal_transition_out_dir, recursive=TRUE)
        setwd(legal_transition_out_dir)
        Spectre::write.files(acc$legal_transitions, i)
        
        # write out the legal and non legal transitions
        illegal_transition_out_dir <- paste0(output.dir, "/fine/illegal_trans_", option)
        dir.create(illegal_transition_out_dir, recursive=TRUE)
        setwd(illegal_transition_out_dir)
        Spectre::write.files(acc$illegal_transitions, i)
        
        all.scores <- rbind(all.scores, c(parameter = i, ari, tracking = acc$tracking_accuracy))
      } 
   }
  }
  setwd(output.dir)
  Spectre::write.files(all.scores, paste0('all_scores_tracksom_option', option))
  
  option <- option + 1
}
