#' Calculate tracking accuracy
#' 
#' @description Given a clustered time series data and a set of transition rule, calculate the number of transitions that are legal.
#' The metric is adapted from tracking accuracy used in ChronoClust publication.
#' @seealso https://github.com/ghar1821/Chronoclust
#' 
#' @param dat Data table storing the clustered gated dataset
#' @param timepoints Vector containing the time points in order
#' @param timepoint_col Character Column denoting time points in dat
#' @param cluster_lineage_col Character Column denoting cluster id produced by tracking by lineage determination in dat
#' @param cluster_assoc_col Character Column denoting cluster id produced by tracking by association in dat
#' @param population_col Character Column denoting cell phenotype (population) in dat. This represents the ground truth of the cells
#' @param rule Data frame or Data table containing the transition rule. must contain 2 columns (from and to)
#' 
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
#' rule = data.table(
#'   from=c("T cell", "B cell", 'Eosinophil'),
#'   to=c("T cell", "B cell", 'Eosinophil')
#' )
#' 
#' tracking_acc <- calculate_tracking_accuracy(dat = dataset,
#'                                             timepoints = c(1:2),
#'                                             timepoint_col = "timepoint",
#'                                             cluster_lineage_col = "lineage",
#'                                             cluster_assoc_col = "association",
#'                                             population_col = "population",
#'                                             rule = rule)
#'                                             
#' @author Givanna Putri
#' @export

calculate_tracking_accuracy <- function(dat, timepoints, timepoint_col, 
                                        cluster_lineage_col, population_col,
                                        rule) {
  require(data.table)
  
  # dat = cell.dat
  # rule = rule
  # timepoints = unique(cell.dat$FileName)
  # timepoint_col = 'FileName'
  # cluster_lineage_col = 'TrackSOM_metacluster_lineage_tracking'
  # cluster_assoc_col = 'TrackSOM_metacluster_proximity_tracking'
  # population_col = 'PopName'
  
  # TODO remove me.
  # dat = solution.list 
  # timepoints = c(1:5) 
  # timepoint_col = 'FileNo' 
  # cluster_lineage_col = 'TrackSOM_metacluster_lineage_tracking'
  # population_col = 'PopName'
  # rule=transitions
  
  rule_df <- as.data.frame(rule)
  
  # compute cluster mapping
  cluster_mapping <- map_clusters(dat = dat, 
                                  timepoint_col = timepoint_col, 
                                  timepoints = timepoints,
                                  cluster_col = cluster_lineage_col,
                                  population_col = population_col)
  
  # get cluster transition
  cluster_transitions <- get_cluster_transitions(dat, timepoints, timepoint_col, cluster_lineage_col, cluster_assoc_col)
  cluster_transitions <- data.table(cluster_transitions)
  
  # print(cluster_transitions)
  
  num_transitions <- nrow(cluster_transitions)
  n_legal_transitions <- 0
  
  pop_col <- c(population_col)
  
  legal_transitions <- list()
  illegal_transitions <- list()
  
  legal_row <- 1
  illegal_row <- 1
  
  for (timepoint_idx in 2:length(timepoints)) {
    curr_timepoint <- timepoints[timepoint_idx]
    prev_timepoint <- timepoints[timepoint_idx-1]
    
    # Don't really need to filter with both from and to, but for safety.
    transitions <- cluster_transitions[(cluster_transitions$from_timepoint == prev_timepoint &
                                         cluster_transitions$to_timepoint == curr_timepoint),]
    
    mappings <- cluster_mapping[cluster_mapping[[timepoint_col]] %in% c(prev_timepoint, curr_timepoint)]
    
    for (i in 1:nrow(transitions)) {
      transition <- transitions[i,]
      from_population <- mappings[(mappings[[cluster_lineage_col]] == transition$from &
                                    mappings[[timepoint_col]] == transition$from_timepoint),][[pop_col]]
      to_population <- mappings[(mappings[[cluster_lineage_col]] == transition$to &
                                    mappings[[timepoint_col]] == transition$to_timepoint),][[pop_col]]
      
      transition_population <- data.frame(from=c(as.character(from_population)), 
                                          to=c(as.character(to_population)))
      is_transition_legal <- nrow(merge(rule_df, transition_population)) > 0
      
      transition_row <- data.table(from_tp=timepoints[timepoint_idx-1],
                                   from_cluster=transition$from,
                                   from_population=from_population,
                                   to_tp=timepoints[timepoint_idx],
                                   to_cluster=transition$to,
                                   to_population=to_population)

      if(is_transition_legal) {
        n_legal_transitions <- n_legal_transitions + 1
        legal_transitions[[legal_row]] <- transition_row
        legal_row <- legal_row + 1
      } 
      else {
        illegal_transitions[[illegal_row]] <- transition_row
        illegal_row <- illegal_row + 1
      }
    }
  }
  
  # print(n_legal_transitions)
  # print(num_transitions)
  
  legal_transitions <- rbindlist(legal_transitions)
  illegal_transitions <- rbindlist(illegal_transitions)
  
  to_return <- list()
  to_return[['legal_transitions']] <- legal_transitions
  to_return[['illegal_transitions']] <- illegal_transitions
  to_return[['tracking_accuracy']] <- n_legal_transitions/num_transitions
  
  return(to_return)
}


#' Get cluster transitions
#' 
#' @description Given a clustered time series data, parse the cluster ID and work out the transitions
#' 
#' @param dat Data table storing the clustered gated dataset
#' @param timepoints Vector containing the time points in order
#' @param timepoint_col Character Column denoting time points in dat
#' @param cluster_lineage_col Character Column denoting cluster id produced by tracking by lineage determination in dat
#' @param cluster_assoc_col Character Column denoting cluster id produced by tracking by association in dat
#' 
#' @author Givanna Putri
#' @export

get_cluster_transitions <- function(dat, timepoints, timepoint_col, cluster_lineage_col, cluster_assoc_col) {
  
  require(gtools)
  require(stringr)
  require(data.table)
  
  # Remove cells not assigned clusters
  dat_bk <- dat[!is.na(dat[[cluster_lineage_col]]),]
  
  edge_df <- data.frame(from=character(),
                        to=character())

  #### Find all the transitions ####
  for (tp_idx in c(2:length(timepoints))) {
    prev_tp_dat <- dat_bk[dat_bk[[timepoint_col]] == timepoints[tp_idx-1], ]
    prev_tp_clust <- mixedsort(unique(prev_tp_dat[[cluster_lineage_col]]))
    
    complex_clusters <- prev_tp_clust[lapply(prev_tp_clust, get_idx_roundclsbracket) > -1]
    simple_n_pipe_clusters <- setdiff(prev_tp_clust, complex_clusters)
    pipe_clusters <- simple_n_pipe_clusters[sapply(simple_n_pipe_clusters, get_idx_pipe) > -1]
    simple_clusters <- setdiff(simple_n_pipe_clusters, pipe_clusters)
    
    # any of the above could have empty character in it which we don't want.
    # this can happen when no merging is allowed. Complex_clusters is empty
    complex_clusters <- complex_clusters[!complex_clusters %in% ""]
    pipe_clusters <- pipe_clusters[!pipe_clusters %in% ""]
    simple_clusters <- simple_clusters[!simple_clusters %in% ""]
    
    # Order the vector based on the length of each element
    simple_clusters <- simple_clusters[order(nchar(simple_clusters), simple_clusters, decreasing = TRUE)]
    pipe_clusters <- pipe_clusters[order(nchar(pipe_clusters), pipe_clusters, decreasing = TRUE)]
    complex_clusters <- complex_clusters[order(nchar(complex_clusters), complex_clusters, decreasing = TRUE)]
    
    curr_tp_dat <- dat_bk[dat_bk[[timepoint_col]] == timepoints[tp_idx], ]
    curr_tp_clust <- mixedsort(unique(curr_tp_dat[[cluster_lineage_col]]))
    
    # this filter out clusters that only exist in current time point
    # these are new clusters and have no predecessors
    curr_tp_clust_existing <- sapply(as.vector(curr_tp_clust), function(cl) {
      if (grepl(",", cl, fixed = TRUE) || grepl("|", cl, fixed = TRUE) || cl %in% simple_clusters) {
        return(cl)
      } 
    })
    curr_tp_clust_existing <- unlist(curr_tp_clust_existing)
    
    ## add edge between new cluster in this time point and the association
    ## isolate the new clusters, then obtain all cells belonging to them_
    ## get unique combination of cluster lineage and association
    ## create edge between them
    # TODO: not required. Remove when all is ok.
    # curr_tp_clust_new <- setdiff(curr_tp_clust, curr_tp_clust_existing)
    # dat_curr_tp_clust_new <- curr_tp_dat[curr_tp_dat[[cluster_lineage_col]] %in% curr_tp_clust_new]
    # 
    # cols_cluster_ids <- c(cluster_lineage_col, cluster_assoc_col)
    # unique_dat_curr_tp_clust_new <- unique(dat_curr_tp_clust_new[,cols_cluster_ids, with=FALSE])
    # 
    # if (nrow(unique_dat_curr_tp_clust_new) > 0) {
    #   for (r in 1:nrow(unique_dat_curr_tp_clust_new)) {
    #     row <- unique_dat_curr_tp_clust_new[r,]
    #     assoc_id <- unlist(str_split(row[[cluster_assoc_col]], "&"))
    #     for (id in assoc_id) {
    #       trimmed_id <- str_trim(id)
    #       df <- data.frame(paste0(tp_idx-1,'_', trimmed_id), 
    #                        paste0(tp_idx, '_', row[[cluster_lineage_col]]))
    #       names(df) <- c("from","to")
    #       edge_df <- rbind(edge_df, df)
    #     }
    #   }
    # }
    
    for (cl in curr_tp_clust_existing) {
      clean_cl <- cl
      ## Match the complex clusters first
      for (cls in complex_clusters) {
        cls_found <- grepl(cls, clean_cl, fixed = TRUE)
        if (cls_found) {
          df <- data.frame(paste0(tp_idx-1,'_', cls), paste0(tp_idx, '_', cl))
          names(df) <- c("from","to")
          edge_df <- rbind(edge_df, df)
          clean_cl <- gsub(cls, "", clean_cl, fixed = TRUE)
        }
      }
      
      ## Match the split simple clusters first
      for (cls in pipe_clusters) {
        cls_found <- grepl(cls, clean_cl, fixed = TRUE)
        if (cls_found) {
          
          # To prevent matching something like RR|1 to R|1, we replace R|1 from RR|1 and if there is alphabet, then it's not the right match
          clean_cl_without_cls <- gsub(cls, "", clean_cl, fixed = TRUE)
          character_in_clean_cl_without_cls <- grepl("^[A-Za-z]+$", clean_cl_without_cls, perl = T)
          
          if (!character_in_clean_cl_without_cls) {
            df <- data.frame(paste0(tp_idx-1,'_', cls), paste0(tp_idx, '_', cl))
            names(df) <- c("from","to")
            edge_df <- rbind(edge_df, df)
            clean_cl <- gsub(cls, "", clean_cl, fixed = TRUE)
          }
          
          df <- data.frame(paste0(tp_idx-1,'_', cls), paste0(tp_idx, '_', cl))
          names(df) <- c("from","to")
          edge_df <- rbind(edge_df, df)
          clean_cl <- gsub(cls, "", clean_cl, fixed = TRUE)
        }
      }
      
      ## Then simple clusters
      for (cls in simple_clusters) {
        cls_found <- grepl(cls, clean_cl, fixed = TRUE)
        if (cls_found) {
          df <- data.frame(paste0(tp_idx-1,'_',cls), paste0(tp_idx,'_',cl))
          names(df) <- c("from","to")
          edge_df <- rbind(edge_df, df)
          clean_cl <- gsub(cls, "", clean_cl, fixed = TRUE)
        }
      }
    }
  }
  
  # split up the to and from time point into different column
  edge_df_to_return <- data.frame(
    from_timepoint=timepoints[as.numeric(str_split_fixed(edge_df$from,"_", n=2)[,1])],
    from = str_split_fixed(edge_df$from,"_", n=2)[,2],
    to_timepoint=timepoints[as.numeric(str_split_fixed(edge_df$to,"_", n=2)[,1])],
    to = str_split_fixed(edge_df$to,"_", n=2)[,2],
    stringsAsFactors = FALSE)
  
  return(edge_df_to_return)
  
}

#' Get round bracket
#' 
#' @description Helper method to detect round brackets in cluster ID
#' 
#' @param cl Character cluster ID
#' 
#' @author Givanna Putri
#' @export
get_idx_roundclsbracket <- function(cl) {
  tail(gregexpr("\\)", cl)[[1]], n=1)
}

#' Get pipe character
#' 
#' @description Helper method to detect pipe character in cluster ID
#' 
#' @param cl Character cluster ID
#' 
#' @author Givanna Putri
#' @export
get_idx_pipe <- function(cl) {
  tail(gregexpr("\\|", cl)[[1]], n=1)
}