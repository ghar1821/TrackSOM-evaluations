library(Spectre)
library(data.table)

# Plot Figure 5
data_dir <- '~/Documents/phd/tracksom/inspection_wnv_bm/10/'
work_dir <- '~/Documents/phd/tracksom/inspection_wnv_bm/10/'

setwd(data_dir)

dat = read.files()

dat <- lapply(c(0:7), function(i) {
  fread(paste0("Result_FileName_WNV_D", i, ".csv"))
})
dat <- rbindlist(dat)

gated_dat = read.files(
  "~/Documents/GitHub/FlowSOM-tracking/experiment_materials/wnv_bm_dataset/dataset/gated/"
)
gated_dat = do.merge.files(gated_dat)

timepoints = c('Pre-infection', sapply(c(1:7), function(d)
  paste("WNV Day", d)))
names(timepoints) = c(1:8)

markers.cols = names(dat)[c(1:14)]

setwd(work_dir)

# attach the population
dat$population = gated_dat$PopName

# remove cells not gated. Not interested
dat = dat[dat$population != 'None', ]

dat[, timepoint := timepoints[FileNo]]

# Draw network plot and heatmap ====
## Source files

setwd(data_dir)
mapping = fread("mapping/10.csv")
mapping[, timepoint := timepoints[FileNo]]
mapping = mapping[, c("TrackSOM_metacluster_lineage_tracking",
                      "PopName",
                      "timepoint")]
setnames(mapping, 'PopName', 'population')



setwd(work_dir)
dir.create('timeseries_inferno')
setwd("timeseries_inferno/")

DrawTimeseriesHeatmap(
  dat = dat,
  timepoints = timepoints,
  timepoint.col = 'timepoint',
  cluster.col = 'TrackSOM_metacluster_lineage_tracking',
  # markers = markers.cols,
  markers = "BV510.SCA.1.",
  file.format = 'png',
  font.size = 18,
  legend.position = 'none',
  colours = 'inferno',
  plot.width = 8,
  plot.height = 10,
  min.node.size = 5,
  max.node.size = 10,
  load.temp.data = FALSE,
  arrow.length = 5,
  line.width = 0.5
)

rule = fread(
  "~/Documents/GitHub/FlowSOM-tracking/experiment_materials/wnv_bm_dataset/dataset/transition_rules.csv"
)

setwd(work_dir)
dir.create("timeseries_population/")
setwd("timeseries_population/")

source(
  "~/Documents/GitHub/TrackSOM-evaluations/tracksom_paper/visualisation_scripts/custom_timeseries_heatmap_WNV_BM.R"
)
colours <- fread("custom_colors.csv")
colours_vec <- colours$colour
names(colours_vec) <- colours$population
make.custom.timeseries.heatmap.WNV.BM(
  dat = dat,
  timepoints = timepoints,
  timepoint.col = 'timepoint',
  cluster.col = 'TrackSOM_metacluster_lineage_tracking',
  file.format = 'png',
  font.size = 18,
  legend.position = 'none',
  plot.width = 8,
  plot.height = 10,
  min.node.size = 5,
  max.node.size = 10,
  load.temp.data = FALSE,
  arrow.length = 5,
  line.width = 0.5,
  mapping = mapping,
  population.col = 'population',
  rule = rule,
  custom.colours = TRUE,
  colours = colours_vec
)

for (col in c("spectral", 'viridis', 'inferno')) {
  setwd(work_dir)
  dir.create(paste0("network_plot_", col))
  setwd(paste0("network_plot_", col))
  DrawNetworkPlot(
    dat = dat,
    timepoints = timepoints,
    timepoint.col = 'timepoint',
    cluster.col = 'TrackSOM_metacluster_lineage_tracking',
    marker.cols = "BV510.SCA.1.",
    arrow.head.gap = 3,
    arrow.length = 2,
    img.height = 10,
    img.width = 10,
    min.node.size = 5,
    max.node.size = 10,
    load.temp.data = FALSE,
    mapping = mapping,
    population.col = 'population',
    standard.colours = col,
    graph.layout = 'gem',
    calculation.type = 'mean',
    legend.position = 'bottom',
    line.width = 0.5
  )
}

for (col in c("spectral", 'bupu', 'viridis', 'inferno')) {
  setwd(work_dir)
  dir.create(paste0("network_plot_", col))
  setwd(paste0("network_plot_", col))
  make.network.plot(
    dat = dat,
    timepoints = timepoints,
    timepoint.col = 'timepoint',
    cluster.col = 'TrackSOM_metacluster_lineage_tracking',
    marker.cols = markers.cols,
    arrow.head.gap = 7,
    arrow.length = 4,
    img.height = 30,
    img.width = 30,
    min.node.size = 10,
    max.node.size = 25,
    load.temp.data = TRUE,
    mapping = mapping,
    population.col = 'population',
    standard.colours = col,
    graph.layout = 'gem',
    calculation.type = 'mean'
  )
}
