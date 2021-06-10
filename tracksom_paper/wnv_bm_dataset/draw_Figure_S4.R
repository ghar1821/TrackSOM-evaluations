# investigate chronoclust ====
work_dir = "/Users/givanna/Documents/phd/tracksom/wnv_bm_timeseries/eval_20201021/chronoclust/"
data_dir = "/Users/givanna/Documents/phd/tracksom/wnv_bm_timeseries/eval_20201021/chronoclust/result_export/"
cc_data = Spectre::read.files(data_dir)

count_per_result = list()
for (i in names(cc_data)) {
    dat = cc_data[[i]]
    count = dat[, .(count=uniqueN(tracking_by_lineage)), by=timepoint]
    count$parameter = rep(i, nrow(count))
    count_per_result[[i]] = count
}
count_per_result = rbindlist(count_per_result)
timepoints = c('Pre-infection', sapply(c(1:7), function(d) paste("WNV Day", d)))
names(timepoints) = sapply(c(0:7), function(x) as.character(x))
setnames(count_per_result, "timepoint", "original_cc_timepoint")
count_per_result[,timepoint:= timepoints[as.character(original_cc_timepoint)]]

counts = count_per_result[, .(median=median(count), avg=mean(count), min=min(count), max=max(count)), by=timepoint]
fwrite(counts, "cluster_count_stats.csv")

font.size = 10
x = ggplot(count_per_result, aes(x = timepoint, y = count, fill = timepoint)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    scale_fill_brewer(palette = "Dark2") + 
    # geom_jitter(shape=16, position=position_jitter(0.2)) +
    theme_bw() +
    scale_y_continuous(limits = c(0,4500), breaks = seq(0,4500, by = 500)) +
    # # scale_x_continuous(limits = c(0.35, 1), breaks = seq(0.35, 1, by = 0.1)) +
    xlab("Time-point") +
    ylab("Number of clusters") +
    labs(fill = 'Time-point') +
    theme(text = element_text(size=font.size, color = 'black'),
          legend.position = 'none',
          axis.text.x = element_text(colour = "black", size = font.size, angle = 90, vjust = 0.3, hjust = 1), 
          axis.text.y = element_text(colour = "black", size = font.size), 
          axis.title.x = element_text(colour = "black", size = font.size + 2),
          axis.title.y = element_text(colour = "black", size = font.size + 2),
          title = element_text(colour = "black", size = font.size + 2),
          rect = element_rect(size = 1),
          panel.grid.minor = element_line(size = 1),
          panel.grid.major = element_line(size = 1)) +
    ggtitle("Number of clusters produced by ChronoClust per time-point")
setwd(work_dir)
ggsave(filename = "number_of_clusters.pdf",
       plot = x,
       dpi = 1200,
       width = 149,
       height = 210,
       units = 'mm')
count_sensible = count_per_result[count_per_result$count >= 10 & count_per_result$count <= 40,]
