# Draw plots and run PRCC

library(ggplot2)
library(data.table)
library(gtools)
library(Spectre)

work_dir = '~/Documents/phd/tracksom/inspection_wnv_bm/20210519/result_qualities/'
setwd(work_dir)

quality_dat_list = read.files()
quality_dat = Spectre::do.merge.files(quality_dat_list)

# read lhc samples and find duplicates
lhc_samples = read.files("~/Documents/GitHub/FlowSOM-tracking/experiment_materials/wnv_bm_dataset/lhc_samples")
lhc_tracksom = lhc_samples[['TrackSOM_lhc']]

# round the numbers so we know which is duplicated
cols = names(lhc_tracksom)[1:10]
lhc_tracksom[,(cols) := round(.SD), .SDcols=cols]

# find duplicates based on just max meta and grid size (for aa and pi mode)
lhc_tracksom$dup_for_AA_PI = duplicated(lhc_tracksom, by = names(lhc_tracksom)[c(1,10)])
lhc_tracksom$dup_for_PV = duplicated(lhc_tracksom, by = names(lhc_tracksom)[c(2:10)])

dup_aa_pi = lhc_tracksom[lhc_tracksom$dup_for_AA_PI == TRUE,]$parameter
dup_pv = lhc_tracksom[lhc_tracksom$dup_for_PV == TRUE,]$parameter # this is nothing! good!

# remove duplicates in quality dat
# 1 and 2 is AA and PI
quality_dat_noDup = quality_dat[!(quality_dat$parameter %in% dup_aa_pi &
                                  quality_dat$FileName %in% c("with_merge_aa", "with_merge_pi", "no_merge_aa", "no_merge_pi"))]

# count number of solutions per operation
quality_dat_noDup[, .(count = .N), by = FileName]

# map the naming properly
modes = c("ChronoClust",
          "No Merging Autonomous Adaptive", "No Merging Prescribed Invariant", "No Merging Prescribed Variant",
          "Merging Autonomous Adaptive", "Merging Prescribed Invariant", "Merging Prescribed Variant")
names(modes) = unique(quality_dat$FileName)
# check this is ok
modes
quality_dat_noDup[, mode:=modes[FileName]]

# plot ecdf and do ks test ====

min_ari <- round(min(quality_dat_noDup$ari_mean), 1) - 0.1
max_ari <- round(max(quality_dat_noDup$ari_mean), 1) + 0.1
min_tracking <- round(min(quality_dat_noDup$tracking), 1) - 0.1
max_tracking <- round(max(quality_dat_noDup$tracking), 1)

# Draw ecdf tracking accuracy ----
df.plt <- ggplot(quality_dat_noDup, aes(tracking, colour = mode)) + 
  stat_ecdf(geom = "step") +
  labs(x=NULL, y=NULL) +
  scale_color_manual(breaks = modes,
                     values = c("#cccccc", "#4b9c7a", "#ca6628", "#7470af", "#7ec0a6", "#ed926b", "#929fc7")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(min_tracking, max_tracking), breaks = seq(min_tracking, max_tracking, by = 0.1)) +
  theme_bw() +
  theme(text = element_text(size=20),
        legend.position = 'none',
        rect = element_rect(size = 1),
        panel.grid.minor = element_line(size = 1), 
        panel.grid.major = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = grid::unit(c(0,0,0,0), "mm"))
ggsave(filename = "ecdf_tracking.png",
       plot = df.plt,
       dpi = 1200,
       width = 150,
       height = 100,
       units = 'mm')
df.plt <- df.plt + theme(legend.position = 'right')
ggsave(filename='ecdf_tracking_legend.png', width=300, height=100, units='mm')

# draw ecdf ARI ----
df.plt <- ggplot(quality_dat_noDup, aes(ari_mean, colour = mode)) + 
  stat_ecdf(geom = "step") +
  labs(x=NULL, y=NULL) +
  scale_color_manual(breaks = modes,
                     values = c("#cccccc", "#4b9c7a", "#ca6628", "#7470af", "#7ec0a6", "#ed926b", "#929fc7")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_x_continuous(limits = c(min_ari, max_ari), breaks = seq(min_ari, max_ari, by = 0.1)) +
  theme_bw() +
  theme(text = element_text(size=20),
        legend.position = 'none',
        rect = element_rect(size = 1),
        panel.grid.minor = element_line(size = 1), 
        panel.grid.major = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = grid::unit(c(0,0,0,0), "mm"))
ggsave(filename = "ecdf_ari.png",
       plot = df.plt,
       dpi = 1200,
       width = 150,
       height = 100,
       units = 'mm')
df.plt <- df.plt + theme(legend.position = 'right')
ggsave(filename='ecdf_ari_legend.png', width=300, height=100, units='mm')

## plot scatter ====
b <- ggplot(quality_dat_noDup, aes_string(x = "ari_mean", y = "tracking")) +
  geom_point(aes_string(color = "mode", shape = 'mode'), size=1) +
  theme_bw() +
  scale_color_manual(breaks = modes,
                     values = c("#cccccc", "#4b9c7a", "#ca6628", "#7470af", "#7ec0a6", "#ed926b", "#929fc7")) +
  scale_shape_manual(breaks = modes,
                     values = c(15, 17, 17, 17, 19, 19, 19)) + 
  scale_y_continuous(limits = c(min_tracking, max_tracking), breaks = seq(min_tracking, max_tracking, by = 0.1)) +
  scale_x_continuous(limits = c(min_ari, max_ari), breaks = seq(min_ari, max_ari, by = 0.1)) +
  xlab("Mean ARI") +
  ylab("Tracking Accuracy") +
  theme(text = element_text(size=30),
        legend.position = 'right',
        rect = element_rect(size = 1),
        panel.grid.minor = element_line(size = 1),
        panel.grid.major = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = grid::unit(c(0,0,0,0), "mm"))
ggsave(filename = "scatter_ari_tracking_legend.png",
       plot = b,
       dpi = 1200,
       width = 600,
       height = 100,
       units = 'mm')
b <- b + theme(legend.position = 'none')
ggsave(filename = "scatter_ari_tracking.png",
       plot = b,
       dpi = 1200,
       width = 200,
       height = 100,
       units = 'mm')

# Compute KS test ----

fileno_combo = gtools::combinations(7, 2, unique(quality_dat_noDup$mode))

ks_res_tracking = list()
for (row in c(1:nrow(fileno_combo))) {
  combo_1 = fileno_combo[row, 1]
  combo_2 = fileno_combo[row, 2]
  
  tracking_combo_1 = quality_dat_noDup[quality_dat_noDup$mode == combo_1]$tracking
  tracking_combo_2 = quality_dat_noDup[quality_dat_noDup$mode == combo_2]$tracking
  
  ks_val = ks.test(tracking_combo_1, tracking_combo_2)
  
  ks_res_tracking[[paste(combo_1, combo_2, sep=',')]] = ks_val
}

ks_res_ari = list()
for (row in c(1:nrow(fileno_combo))) {
  combo_1 = fileno_combo[row, 1]
  combo_2 = fileno_combo[row, 2]
  
  ari_combo_1 = quality_dat_noDup[quality_dat_noDup$mode == combo_1]$ari_mean
  ari_combo_2 = quality_dat_noDup[quality_dat_noDup$mode == combo_2]$ari_mean
  
  ks_val = ks.test(ari_combo_1, ari_combo_2)
  
  ks_res_ari[[paste(combo_1, combo_2, sep=',')]] = ks_val
}

ks_res = data.table(
  combo=c(names(ks_res_ari), names(ks_res_tracking)),
  metric=c(rep("mean_ari", 3), rep("tracking", 3)),
  d_val=c(sapply(ks_res_ari, function(r) r$statistic),
          sapply(ks_res_tracking, function(r) r$statistic)),
  p_val=c(sapply(ks_res_ari, function(r) r$p.value),
          sapply(ks_res_tracking, function(r) r$p.value))
)
ks_res[, stat_sig:=p_val < 0.005]
ks_res[, val := ifelse(p_val < 0.005, paste0(format(round(d_val, 2), nsmall=2), '*'), format(round(d_val, 2), nsmall=2))]
ks_res[, c("mode_1", "mode_2") := tstrsplit(combo, ",", fixed=TRUE)]
keep_cols = c('metric', 'mode_1', 'mode_2', 'val')
ks_res = ks_res[, ..keep_cols]
ks_res$mode_1 = factor(ks_res$mode_1, levels = modes)
ks_res$mode_2 = factor(ks_res$mode_2, levels = modes)

# dimnames: column first, then row
ks_res_mat_ari = matrix(nrow=7, ncol=7, dimnames = list(modes, modes))
ks_res_mat_tracking = matrix(nrow=7, ncol=7, dimnames = list(modes, modes))

# assume acessing, row then column.
for (i in c(1:nrow(ks_res))) {
  res_row = ks_res[i,]
  
  if (res_row$metric == 'mean_ari') {
    row_idx = match(res_row$mode_1, modes)
    col_idx = match(res_row$mode_2, modes)
    ks_res_mat_ari[col_idx, row_idx] = res_row$val
  } else {
    col_idx = match(res_row$mode_1, modes)
    row_idx = match(res_row$mode_2, modes)
    ks_res_mat_tracking[col_idx, row_idx] = res_row$val
  }
}
# have to move the merging ones to the bottom
# view the matrix, and you will know why!
for (row in c(2:4)) {
  for (col in c(5:7)) {
    ks_res_mat_ari[col, row] = ks_res_mat_ari[row, col]
    ks_res_mat_ari[row, col] = NA
  }
}

for (row in c(5:7)) {
  for (col in c(2:4)) {
    ks_res_mat_tracking[col, row] = ks_res_mat_tracking[row, col]
    ks_res_mat_tracking[row, col] = NA
  }
}

combined_ks_res_mat = matrix(nrow=7, ncol=7, dimnames = list(modes, modes))
combined_ks_res_mat[upper.tri(combined_ks_res_mat)] = ks_res_mat_tracking[upper.tri(ks_res_mat_tracking)]
combined_ks_res_mat[lower.tri(combined_ks_res_mat)] = ks_res_mat_ari[lower.tri(ks_res_mat_ari)]

write.table(combined_ks_res_mat, file='ks_matrix.csv', sep = ",")
fwrite(ks_res, "ks_test.csv")

# PRCC ====

# remove rounding of lhc samples

quality_dat_tracksom = quality_dat_noDup[quality_dat_noDup$mode != 'ChronoClust']
quality_dat_tracksom[, .(count = .N), by = mode]
quality_dat_tracksom = merge(quality_dat_tracksom, lhc_tracksom, by = 'parameter', suffixes = c("_from_data", "_from_param"))
cols = names(quality_dat_tracksom)

mode_op = c("Autonomous Adaptive", "Prescribed Invariant", "Prescribed Variant")
names(mode_op) = c("AA","PI","PV")
merge_op = c("Merging", "No Merging")
names(merge_op) = c("M", "NM")

library(spartan)

# calculate PRCC score
setwd(work_dir)
dir.create("prcc_20210520")
setwd("prcc_20210520")

prcc_scores_list = list()
for (i in c(1:length(mode_op))) {
  for (j in c(1:length(merge_op))) {
    mo = mode_op[i]
    me = merge_op[j]
    quality_dat_tracksom_subset = quality_dat_tracksom[quality_dat_tracksom$mode == paste(me, mo)]
    corcoefffilename = gsub(" ", "", paste("LHC_corCoeffs_",names(me),"_",names(mo), ".csv"))

    if (i %in% c(1,2)) {
      params = cols[c(15,24)]
    } else {
      params = cols[c(16:24)]
    }

    fname = paste0(names(me),"_",names(mo), ".csv")
    fwrite(quality_dat_tracksom_subset, fname)

    prcc_val <- lhc_generatePRCoEffs(FILEPATH = getwd(),
                         PARAMETERS = params,
                         MEASURES = c("ari_mean","tracking"),
                         LHCSUMMARYFILENAME = fname,
                         CORCOEFFSOUTPUTFILE = corcoefffilename,
                         cor_calc_method = 's')
    score = fread(corcoefffilename)
    score$mode = paste(names(me), names(mo), sep = '_')
    score$clustering_mode = paste0(mo, " (", names(mo), ")")
    score$merging_mode = paste0(me, " (", names(me), ")")
    setnames(score, 'V1', 'parameter')
    prcc_scores_list[[paste(names(me), names(mo))]] = score
  }
}

prcc_scores_list <- rbindlist(prcc_scores_list)

prcc_scores_statsig_ari <- prcc_scores_list[prcc_scores_list$ari_mean_PValue < 0.005, c("parameter", "mode", 
                                                                                        "clustering_mode", "merging_mode",
                                                                                        "ari_mean_Estimate", "ari_mean_PValue")]
prcc_scores_statsig_tracking <- prcc_scores_list[prcc_scores_list$tracking_PValue < 0.005,
                                                 c("parameter", "mode", 
                                                   "clustering_mode", "merging_mode",
                                                   "tracking_Estimate", "tracking_PValue")]

fwrite(prcc_scores_statsig_ari, "sig_ari.csv")
fwrite(prcc_scores_statsig_tracking, "sig_tracking.csv")

library(ggplot2)
library(scales)

modes_to_plot <- c("M_AA", "NM_AA", "M_PI", "NM_PI", "M_PV",  "NM_PV")
names(modes_to_plot) <- c("Merging Autonomous Adaptive",
                          "No Merging Autonomous Adaptive",
                          "Merging Prescribed Invariant",
                          "No Merging Prescribed Invariant",
                          "Merging Prescribed Variant",
                          "No Merging Prescribed Variant")

for (mode_to_plot in modes_to_plot) {
  prcc_to_plot <- fread(paste0(mode_to_plot, ".csv"))
  params_to_plot_ari <- prcc_scores_statsig_ari[prcc_scores_statsig_ari$mode == mode_to_plot,]$parameter

  for (param_to_plot_ari in params_to_plot_ari) {
    prcc_to_plot_ari <- prcc_to_plot[, c("ari_mean", param_to_plot_ari), with=FALSE]

    plt <- ggplot(prcc_to_plot_ari, aes_string(y="ari_mean", x=param_to_plot_ari)) +
      geom_point() +
      scale_x_continuous(breaks=pretty_breaks(10)) +
      scale_y_continuous(breaks=pretty_breaks(10)) +
      theme(panel.background = element_rect(fill = 'white', colour = 'grey'),
            axis.text = element_text(size = 16),
            axis.title = element_blank(),
            legend.position = 'bottom')
      # ggtitle(paste("PRCC analysis for", names(mode_to_plot))) +
      # labs(x='Parameter values', y="Mean ARI")

    ggsave(filename = paste0("ari_", mode_to_plot, "_",  param_to_plot_ari, ".png"),
           plot = plt,
           width = 5,
           height = 5,
           dpi = 1200)
  }


  params_to_plot_tracking <- prcc_scores_statsig_tracking[prcc_scores_statsig_tracking$mode == mode_to_plot,]$parameter

  for (param_to_plot_tracking in params_to_plot_tracking) {
    prcc_to_plot_tracking <- prcc_to_plot[, c("tracking", param_to_plot_tracking), with=FALSE]
    plt <- ggplot(prcc_to_plot_tracking, aes_string(y="tracking", x=param_to_plot_tracking)) +
      geom_point() +
      scale_x_continuous(breaks=pretty_breaks(10)) +
      scale_y_continuous(breaks=pretty_breaks(10)) +
      theme(panel.background = element_rect(fill = 'white', colour = 'grey'),
            axis.text = element_text(size = 16),
            axis.title = element_blank(),
            legend.position = 'bottom')
      # ggtitle(paste("PRCC analysis for", names(mode_to_plot))) +
      # labs(x='Parameter values', y="Mean ARI")

      ggsave(filename = paste0("tracking_", mode_to_plot, "_",  param_to_plot_tracking, ".png"),
             plot = plt,
             width = 5,
             height = 5,
             dpi = 1200)
  }


}
