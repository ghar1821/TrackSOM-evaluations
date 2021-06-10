#' Map clusters
#' 
#' @description 
#' The function to map clusters to ground truth label (population) based on the highest proportion of cells it captured.
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

map_clusters <- function(dat, timepoint_col, timepoints, cluster_col, population_col) {
  # dat=dataset 
  # timepoint_col="timepoint" 
  # cluster_col="cluster"
  # population_col="population"
  
  require(data.table)

  dat_bk <- data.table(dat)
  
  group_by_col <- c(cluster_col, population_col)
  
  mapping_per_timepoint <- list()
  
  for (t in timepoints) {
    # t = timepoints[1]
    dat_sub <- dat_bk[dat_bk[[timepoint_col]] == t,]
    
    count_per_population_cluster <- dat_sub[, .(count = .N), by=group_by_col]
    mapping <- count_per_population_cluster[ , .SD[which.max(count)], by = eval(cluster_col)]
    mapping[[timepoint_col]] <- rep(t, nrow(mapping))
    mapping_per_timepoint[[t]] <- mapping
  }
  
  mapping_to_return <- rbindlist(mapping_per_timepoint)
  
  return(mapping_to_return)
}