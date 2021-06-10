## Investigate more details on the synthetic data ====
## Miscellaneous scripts and codes ====
library(ggplot2)
library(data.table)
library(gtools)

# count number of cell types per day
gated_dat = read.files(
  "~/Documents/GitHub/FlowSOM-tracking/experiment_materials/synthetic_dataset/dataset/csv_files_withLabel/"
)
gated_dat = do.merge.files(gated_dat)

# plot ecdf and do ks test ====
work_dir = '~/Documents/GitHub/FlowSOM-tracking/experiment_materials/synthetic_dataset/result_qualities/'
setwd(work_dir)

quality_dat_list = Spectre::read.files('.')
quality_dat = Spectre::do.merge.files(quality_dat_list)

# read lhc samples and find duplicates
lhc_samples = fread(
  '~/Documents/GitHub/FlowSOM-tracking/experiment_materials/synthetic_dataset/lhc_samples/LHC_sample.csv'
)
# round the numbers so we know which is duplicated
cols = names(lhc_samples)[1:7]
round_cols = sapply(cols, function(s)
  paste0(s, "_round"))
lhc_samples[, (round_cols) := round(.SD), .SDcols = cols]

# find duplicates based on just max meta and grid size (for aa and pi mode)
lhc_samples$dup_for_AA_PI = duplicated(lhc_samples, by = round_cols[c(1, 7)])
lhc_samples$dup_for_PV = duplicated(lhc_samples, by = round_cols[c(2:7)])

dup_aa_pi = lhc_samples[lhc_samples$dup_for_AA_PI == TRUE, ]$parameter
dup_pv = lhc_samples[lhc_samples$dup_for_PV == TRUE, ]$parameter # this is nothing! good!

# remove duplicates in quality data
# 1 and 2 is AA and PI
quality_dat_noDup = quality_dat[!(quality_dat$parameter %in% dup_aa_pi &
                                    quality_dat$FileNo %in% c(1, 2))]

# count number of solutions per operation
quality_dat_noDup[, .(count = .N), by = FileNo]

# ecdf for tracking accuracy
df.plt <-
  ggplot(quality_dat_noDup, aes(tracking, colour = FileName)) +
  stat_ecdf(geom = "step") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0.84, 1), breaks = seq(0.8, 1, by = 0.02)) +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    legend.position = 'none',
    rect = element_rect(size = 1),
    panel.grid.minor = element_line(size = 1),
    panel.grid.major = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
  )
ggsave(
  filename = "ecdf_tracking.png",
  plot = df.plt,
  dpi = 1200,
  width = 150,
  height = 100,
  units = 'mm'
)

# ecdf for ARI
df.plt <-
  ggplot(quality_dat_noDup, aes(ari_mean, colour = FileName)) +
  stat_ecdf(geom = "step") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0.30, 1), breaks = seq(0, 1, by = 0.1)) +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    legend.position = 'none',
    rect = element_rect(size = 1),
    panel.grid.minor = element_line(size = 1),
    panel.grid.major = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
ggsave(
  filename = "ecdf_ari.png",
  plot = df.plt,
  dpi = 1200,
  width = 150,
  height = 100,
  units = 'mm'
)

# KS test
fileno_combo = gtools::combinations(3, 2, c(1, 2, 3))

ks_res_tracking = list()
for (row in c(1:nrow(fileno_combo))) {
  combo_1 = fileno_combo[row, 1]
  combo_2 = fileno_combo[row, 2]
  
  tracking_combo_1 = quality_dat_noDup[quality_dat_noDup$FileNo == combo_1]$tracking
  tracking_combo_2 = quality_dat_noDup[quality_dat_noDup$FileNo == combo_2]$tracking
  
  ks_val = ks.test(tracking_combo_1, tracking_combo_2)
  
  ks_res_tracking[[paste(combo_1, combo_2)]] = ks_val
}

ks_res_ari = list()
for (row in c(1:nrow(fileno_combo))) {
  combo_1 = fileno_combo[row, 1]
  combo_2 = fileno_combo[row, 2]
  
  ari_combo_1 = quality_dat_noDup[quality_dat_noDup$FileNo == combo_1]$ari_mean
  ari_combo_2 = quality_dat_noDup[quality_dat_noDup$FileNo == combo_2]$ari_mean
  
  ks_val = ks.test(ari_combo_1, ari_combo_2)
  
  ks_res_ari[[paste(combo_1, combo_2)]] = ks_val
}

ks_res = data.table(
  combo = c(names(ks_res_ari), names(ks_res_tracking)),
  metric = c(rep("mean_ari", 3), rep("tracking", 3)),
  d_val = c(
    sapply(ks_res_ari, function(r)
      r$statistic),
    sapply(ks_res_tracking, function(r)
      r$statistic)
  ),
  p_val = c(
    sapply(ks_res_ari, function(r)
      r$p.value),
    sapply(ks_res_tracking, function(r)
      r$p.value)
  )
)
ks_res[, stat_sig := p_val < 0.005]

fwrite(ks_res, "ks_test.csv")

## plot scatter ====
options <- c("Auto", "Invariant", "Variant")
names(options) <- unique(quality_dat_noDup$FileNo)

quality_dat_noDup[, Mode := options[FileNo]]

b <-
  ggplot(quality_dat_noDup, aes_string(x = "ari_mean", y = "tracking")) +
  geom_point(aes_string(color = "Mode"), size = 1) +
  theme_bw() +
  scale_y_continuous(limits = c(0.75, 1),
                     breaks = seq(0.75, 1, by = 0.05)) +
  scale_x_continuous(limits = c(0.35, 1), breaks = seq(0.35, 1, by = 0.1)) +
  theme(
    text = element_text(size = 20),
    legend.position = 'none',
    rect = element_rect(size = 1),
    panel.grid.minor = element_line(size = 1),
    panel.grid.major = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )
ggsave(
  filename = "scatter_ari_tracking.png",
  plot = b,
  dpi = 1200,
  width = 150,
  height = 100,
  units = 'mm'
)

# Draw network and heatmap plot for single result ====
source(
  "~/Documents/GitHub/TrackSOM-evaluations/tracksom_paper/visualisation_scripts/custom_timeseries_heatmap.R"
)

work_dir = '~/Documents/GitHub/FlowSOM-tracking/experiment_materials/synthetic_dataset/raw_result_in_manuscript/'
setwd(work_dir)
setwd("86")
dat = read.files()
dat = do.merge.files(dat)
setnames(dat, 'x', 'dimX')
setnames(dat, 'y', 'dimY')
setnames(dat, 'z', 'dimZ')
dat[, timepoint := paste("Day", FileNo)]
timepoints = sapply(c(1:5), function(d)
  paste("Day", d))

gating.dir <-
  "~/Documents/GitHub/FlowSOM-tracking/experiment_materials/synthetic_dataset/dataset/csv_files_withLabel/"
gating.list <- Spectre::read.files(
  file.loc = gating.dir,
  file.type = '.csv',
  do.embed.file.names = TRUE
)
gating.dat <- Spectre::do.merge.files(gating.list)

rule = fread(
  "~/Documents/GitHub/FlowSOM-tracking/experiment_materials/synthetic_dataset/dataset/transition_rules.csv"
)

setwd(
  "~/Documents/GitHub/TrackSOM-evaluations/tracksom_paper/synthetic_dataset/figure_3/"
)
mapping = fread("cluster_pop_mapping.csv")
mapping[, timepoint := paste("Day", FileNo)]

cluster_order = fread("custom_order.csv", sep = "")$order
custom_colours <- fread("custom_color.csv")
custom_colours_vector <- custom_colours$Colour
names(custom_colours_vector) <- custom_colours$Population

# TODO add the custom colour functionality to the TrackSOM timeseries heatmap.
make.custom.timeseries.heatmap(
  dat = dat,
  timepoints = timepoints,
  timepoint.col = 'timepoint',
  cluster.col = 'TrackSOM_metacluster_lineage_tracking',
  markers = c('dimX', 'dimY', 'dimZ'),
  file.format = 'png',
  font.size = 20,
  legend.position = 'none',
  min.node.size = 5,
  max.node.size = 20,
  load.temp.data = FALSE,
  arrow.length = 5,
  mapping = mapping,
  population.col = 'PopName',
  rule = rule,
  cluster.ordering = cluster_order,
  line.width = 0.5,
  colours = custom_colours_vector
)

setwd(work_dir)
TrackSOM::DrawTimeseriesHeatmap(
  dat = dat,
  timepoints = timepoints,
  timepoint.col = 'timepoint',
  cluster.col = 'TrackSOM_metacluster_lineage_tracking',
  markers = c('dimX', 'dimY', 'dimZ'),
  file.format = 'png',
  font.size = 20,
  legend.position = 'none',
  min.node.size = 5,
  max.node.size = 20,
  colours = 'spectral',
  cluster.ordering = cluster_order,
  line.width = 0.5
)

TrackSOM::DrawNetworkPlot(
  dat = dat,
  timepoints = timepoints,
  timepoint.col = 'timepoint',
  cluster.col = 'TrackSOM_metacluster_lineage_tracking',
  marker.cols = c('dimX', 'dimY', 'dimZ'),
  min.node.size = 5,
  max.node.size = 15,
  arrow.head.gap = 7,
  standard.colours = 'cividis',
  legend.position = 'none',
  file.format = 'png'
)
