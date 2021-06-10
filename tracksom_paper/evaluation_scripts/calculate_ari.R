#' Compute ARI
#' 
#' @description 
#' The function calculate ARI using mclust
#' 
#' @param dat Data table storing the clustered gated dataset
#' @param timepoints Vector containing the time points in order
#' @param timepoint_col Character Column denoting time points in dat
#' @param cluster_col Character Column denoting clusters in dat
#' @param population_col Character Column denoting cell phenotype (population) in dat. This represents the ground truth of the cells
#' 
#' @return aris List containing ARI per time point as well as average
#' 
#' @author Givanna Putri
#' @example 
#' SCA1=sample(seq(100,400), size=15),
#' timepoint=c(
#'   rep(1,5),
#'   rep(2,10)
#' ),
#' population=c(rep('T cell',3),
#'              rep('B cell', 2),
#'              rep('T cell', 5),
#'              rep('B cell', 2),
#'              rep('Eosinophil', 3)
#' ),
#' lineage=c(rep('A',3),
#'           rep('B', 2),
#'           rep('(A,B)', 5),
#'           rep('B|1',2),
#'           rep('C', 3)
#' ),
#' association=c(rep('None',5),
#'               rep('A&B',5),
#'               rep('B', 2),
#'               rep('A',3)
#' )
#' )
#' 
#' 
#' ari <- calculate_ari(dat = dataset,
#' timepoints = c(1:2),
#' timepoint_col = "timepoint",
#' cluster_col = "lineage",
#' population_col = "population")
#' @export

calculate_ari <- function(dat, timepoints, timepoint_col, cluster_col, population_col) {
  require(mclust)
  
  dat_bk <- data.table(dat)
  
  aris <- list()
  for (t in timepoints) {
    dat_sub <- dat_bk[dat_bk[[timepoint_col]] == t,]
    ari <- mclust::adjustedRandIndex(as.vector(dat_sub[[population_col]]),
                                     as.vector(dat_sub[[cluster_col]]))
    aris[[t]] <- ari
  }
  
  average_ari <- mean(unlist(aris))
  aris[['mean']] <- average_ari
  
  return(aris)
}