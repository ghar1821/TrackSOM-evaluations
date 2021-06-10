# custom time series heatmap which colour the nodes by population, 
# shaped the nodes by activated or unactivated
# change line type based on whether the rule make sense.
# just for the paper

library(Spectre)
library(ggplot2)
library(viridis)
library(data.table)
library(stringr)

# setup ----
# dat = dat
# timepoints = timepoints
# timepoint.col = 'timepoint'
# cluster.col = 'FlowSOM_metacluster_lineage_tracking'
# markers = markers.cols
# file.format = 'pdf'
# font.size = 25
# legend.position = 'right'
# plot.width = 30
# plot.height = 30
# min.node.size = 10
# max.node.size = 25
# load.temp.data = FALSE
# arrow.length = 5
# mapping = mapping
# population.col = 'population'
# rule = rule
# cluster.ordering = cluster_order
# map.cluster.id = TRUE

make.custom.timeseries.heatmap.WNV.BM <- function(dat, timepoints, timepoint.col, 
                                           cluster.col, rule,
                                           mapping, population.col,
                                           cluster.ordering = NULL, map.cluster.id = FALSE,
                                           plot.width = 12, plot.height = 15,
                                           font.size = 16, file.format = 'pdf',
                                           min.node.size = 6, max.node.size = 20,
                                           plot.title = NULL, legend.position = 'right',
                                           load.temp.data = FALSE, colours = 'spectral', custom.colours = FALSE,
                                           arrow.length = 2, line.width = 1) {
  
  # custom heatmap where edge are coloured by either transition between
  # different cell types or unactivated -> activated or non-sense
  message("Computing node details")
  
  # So we need a data table with day, cluster id, centroid, proportion column for just the bubble plot
  grp <- c(timepoint.col, cluster.col)
  
  
  cellCount <- dat[, .(count = .N), by = grp]
  cellCount_perDay <- dat[, .(countPerDay = .N), by = timepoint.col]
  cellCount <- merge(cellCount, cellCount_perDay, by=timepoint.col, all.x = TRUE)
  cellCount[, proportion:=(count/countPerDay)]
  
  # append the mapping of cluster and population name.
  dat_chart <- merge(cellCount, mapping, by=grp)
  
  message("Computing edge details")
  
  # Get the edges to link the clusters ====
  # the edges need network plot to work out edges
  edges <- get.transitions.as.edges(dat = dat,
                                    timepoints = timepoints,
                                    timepoint_col = timepoint.col,
                                    cluster_col = cluster.col)
  edges <- data.table(edges)
  # have to convert this to data table with columns: id, timepoint?
  from_tp <- sapply(edges$from, function(fr) {
    tp <- as.numeric(str_split_fixed(fr, "_", 2)[1,1])
    timepoints[tp]
  })
  to_tp <- sapply(edges$to, function(to) {
    tp <- as.numeric(str_split_fixed(to, "_", 2)[1,1])
    timepoints[tp]
  })
  from_cl <- sapply(edges$from, function(fr) {
    str_split_fixed(fr, "_", 2)[1,2]
  })
  to_cl <- sapply(edges$to, function(to) {
    str_split_fixed(to, "_", 2)[1,2]
  })
  
  # colour edges based on legality
  transition_legal = list()
  for (i in c(1:length(to_cl))) {
    from_t <- from_tp[i]
    to_t <- to_tp[i]
    from_clust <- from_cl[i]
    to_clust <- to_cl[i]
    
    from_pop <- mapping[(mapping[[cluster.col]] == from_clust &
                           mapping[[timepoint.col]] == from_t),][[population.col]]
    to_pop <- mapping[(mapping[[cluster.col]] == to_clust &
                         mapping[[timepoint.col]] == to_t),][[population.col]]
    
    rule_exist <- rule[(rule$from == from_pop &
                          rule$to == to_pop),]
    if (nrow(rule_exist) == 1) {
      transition_legal[[i]] = "Plausible"
    } else {
      transition_legal[[i]] = "Not plausible"
    }
  }
  transition_legal <- unlist(transition_legal)
  
  edges_chart <- data.table(timepoint=c(from_tp, to_tp))
  edges_chart$cluster <- c(from_cl,to_cl)
  edges_chart$group <- rep(c(1:length(from_tp)),2)
  edges_chart$transition_valid <- rep(transition_legal, 2)
  
  # order the data based on the timepoints given
  dat_chart[[timepoint.col]] <- factor(dat_chart[[timepoint.col]], levels=timepoints)
  if (!is.null(cluster.ordering)) {
    dat_chart[[cluster.col]] <- factor(dat_chart[[cluster.col]], levels=cluster.ordering)
  } else {
    dat_chart[[cluster.col]] <- factor(dat_chart[[cluster.col]], levels=sort(unique(dat_chart[[cluster.col]])))
  }
  
  # convert ID to numeric
  if (map.cluster.id) {
    if (!is.null(cluster.ordering)) {
      sorted_unique_clust <- cluster.ordering
    } else {
      sorted_unique_clust = sort(unique(dat[[cluster.col]]))
    }
    unique_clust = sapply(c(1: length(sorted_unique_clust)), function(x) as.character(x))
    names(unique_clust) = sorted_unique_clust
    dat_chart$mapped_cluster_id <- sapply(dat_chart[[cluster.col]], function(cl) {
      unique_clust[cl]
    })
    dat_chart$mapped_cluster_id <- factor(dat_chart$mapped_cluster_id, levels=unique_clust)
    edges_chart$mapped_cluster_id <- sapply(edges_chart$cluster, function(cl) {
      unique_clust[cl]
    })
    
    to_plot_node = 'mapped_cluster_id'
    to_plot_edge = 'mapped_cluster_id' 
    
    mapped_dat = data.table(cluster=names(unique_clust),
                            numeric_id=unique_clust)
    setnames(mapped_dat, 'cluster', cluster.col)
    fwrite(mapped_dat, paste0(cluster.col, '_numericMapping.csv'))
  } else {
    to_plot_node = cluster.col
    to_plot_edge = 'cluster' 
  }
  
  # split population to just base
  population_bare <- sapply(dat_chart$population, function(pop) {
    pop_split <- str_split(pop, "_")[[1]]
    popname = paste(head(pop_split, n=length(pop_split)-1), collapse = ' ')
    if (grepl('monocyte', popname)) {
      popname='Monocyte'
    }
    return(popname)
  })
  activation <- sapply(dat_chart$population, function(pop) {
    pop_split <- str_split(pop, "_")[[1]]
    tail(pop_split,1)
  })
  dat_chart$population_bare <- factor(population_bare, levels = mixedsort(unique(population_bare)))
  dat_chart$activation <- activation
  num_populations <- length(unique(dat_chart$population_bare))
  
  message(paste("Drawing timeseries heatmap coloured by", population.col))
  
  # Draw plots ====
  
  xx <-  ggplot(dat_chart, aes_string(x = timepoint.col, y = to_plot_node)) + 
    geom_point(aes(size = proportion, fill = population_bare, shape = factor(activation))) + 
    geom_line(aes_string(x = "timepoint", y = to_plot_edge, group = "group", colour = "transition_valid"), 
              edges_chart, 
              arrow = arrow(length = unit(as.numeric(arrow.length),"mm")), 
              size=line.width) +
    scale_size_continuous(limit = c(min(dat_chart$proportion), max(dat_chart$proportion)),
                          range = c(min.node.size, max.node.size)) +
    labs(x= "", y = "", 
         size = "Proportion of cells", 
         fill = population.col)  + 
    theme(legend.key=element_blank(), 
          axis.text.x = element_text(colour = "black", size = font.size, angle = 90, vjust = 0.3, hjust = 1), 
          axis.text.y = element_text(colour = "black", size = font.size), 
          legend.text = element_text(size = font.size, colour ="black"), 
          legend.title = element_text(size = font.size), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
          legend.position = legend.position,
          plot.title = element_text(color="Black", face="bold", size=font.size, hjust=0),
          plot.margin = grid::unit(c(0,0,0,0), "mm")) +
    scale_shape_manual(values=c(21,24)) +
    # scale_fill_manual(values = c(col.scheme(num_populations))) +
    scale_colour_manual(values=c("red", "black")) +
    guides(colour = guide_legend(override.aes = list(size = 1)),
           shape = guide_legend(override.aes = list(size = min.node.size)),
           fill = guide_legend(override.aes = list(size = min.node.size, shape=21))) +
    scale_y_discrete(limits = rev(levels(factor(dat_chart[[to_plot_node]]))))
  
  if (custom.colours) {
    xx <- xx + scale_fill_manual(values = colours)
  }
  else {
    if (colours == 'spectral') {
      col.scheme <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))
    } else {
      col.scheme <- colorRampPalette(c(viridis_pal(option = colours)(50)))  
    }
    xx <- xx + scale_fill_manual(values = c(col.scheme(num_populations)))
  }
  
  if (map.cluster.id) {
    filename <- paste0("Timeseries_heatmap_mapped_by_", population.col, ".", file.format)
  } else {
    filename <- paste0("Timeseries_heatmap_by_", population.col, ".", file.format)
  }
  
  ggsave(filename = filename,
         plot = xx,
         width = plot.width,
         height = plot.height,
         limitsize = FALSE)
}



