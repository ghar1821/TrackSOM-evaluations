#' Map clusters
#' 
#' @description 
#' The function to map clusters to ground truth label (population) based on their Jaccard similarity.
#' Each cluster is mapped to a population with which it is the most similar to (highest Jaccard similarity score). 
#' 
#' @param dat Data table storing the clustered gated dataset
#' @param timepoint_col Character Column denoting time points in dat
#' @param cluster_col Character Column denoting clusters in dat
#' @param population_col Character Column denoting cell phenotype (population) in dat. This represents the ground truth of the cells
#' 
#' @return mapping Data frame containing mapping
#' 
#' @author Givanna Putri
#' @example 
#' 
#' library(Spectre)
#' setwd("~/Documents/phd/tracksom/sample_covid_data")
#' 
#' dat.list <- Spectre::read.files()
#' dat <- Spectre::do.merge.files(dat.list)
#' 
#' map_clusters(dat, 'FileNo', 'TrackSOM_metacluster_lineage_tracking', 'Population')
#' 
#' @export

compute_cluster_composition <- function(dat, timepoint_col, timepoints, cluster_col, population_col) {
  # dat = dat
  # timepoint_col = timepoint.col
  # timepoints = timepoints
  # cluster_col = 'TrackSOM_metacluster_lineage_tracking'
  # population_col = "population"
  
  require(data.table)
  
  dat_bk <- data.table(dat)
  
  group_by_col <- c(cluster_col, population_col)
  
  mapping_per_timepoint <- list()
  
  for (t in timepoints) {
    # t = timepoints[1]
    dat_sub <- dat_bk[dat_bk[[timepoint_col]] == t,]
    
    count_per_population_cluster <- dat_sub[, .(count = .N), by=group_by_col]
    count_per_cluster <- dat_sub[, .(count = .N), by=c(cluster_col)]
    
    all_count <- merge(count_per_population_cluster, count_per_cluster, by=cluster_col,
                       suffixes = c("_PerPopulationCluster", "PerCluster"))
    all_count[, proportion:=(count_PerPopulationCluster/countPerCluster)]
    all_count[[timepoint_col]] <- rep(t, nrow(all_count))
    mapping_per_timepoint[[t]] <- all_count
  }
  
  mapping_to_return <- rbindlist(mapping_per_timepoint)
  
  return(mapping_to_return)
}