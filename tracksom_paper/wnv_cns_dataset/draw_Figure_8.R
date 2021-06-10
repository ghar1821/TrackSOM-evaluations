library(data.table)
library(Spectre)
library(stringr)
library(scales)
library(TrackSOM)

source("~/Documents/GitHub/TrackSOM-evaluations/tracksom_paper/wnv_cns_dataset/make.autograph.R")

work_dir <- "~/Dropbox (Sydney Uni)/tracksom/wnv_cns/clustered_10x10_pv/"
setwd(work_dir)

dat <- fread("res_10x10_pv.csv")

markers <- names(dat)[c(1:8,10:20)]
timepoints <- c("Mock","WNV-01","WNV-02","WNV-03","WNV-04","WNV-05")

# timeseries heatmap ----
setwd(work_dir)
dir.create('timeseries_inferno')
setwd("timeseries_inferno/")

DrawTimeseriesHeatmap(
  dat = dat,
  timepoints = timepoints,
  timepoint.col = 'Group',
  cluster.col = 'TrackSOM_metacluster_lineage_tracking',
  markers = markers,
  file.format = 'png',
  font.size = 24,
  legend.position = 'right',
  colours = 'inferno',
  plot.width = 12,
  plot.height = 26,
  min.node.size = 6,
  max.node.size = 12,
  load.temp.data = TRUE,
  arrow.length = 5,
  line.width = 1
)

# network plot ----
setwd(work_dir)
dir.create("network_inferno")
setwd("network_inferno")
names(dat) <-
  c(
    "FSC_A",
    "SSC_A",
    "FITC_Ly6C",
    "PerCP_CD45",
    "APC_CD62L",
    "AF700_CD4",
    "APCCy7_CD86",
    "BUV395_CD11b",
    "Live_Dead_Blue",
    "BUV737_B220",
    "BV421_Siglec_F",
    "BV510_IA_IE",
    "BV605_Ly6C",
    "BV650_Ly6G",
    "BV711_CD8a",
    "BV785_CD11c",
    "PE_CD115",
    "PECF594_CD80",
    "PECy5_CD3e",
    "PECy7_NK1_1",
    "FileName",
    "FileNo",
    "Sample",
    "Group",
    "Batch",
    "Cells_per_sample",
    "V6",
    "V7",
    "TrackSOM_cluster",
    "TrackSOM_metacluster",
    "TrackSOM_metacluster_lineage_tracking"
  )
markers <- c(
  "FITC_Ly6C",
  "PerCP_CD45",
  "AF700_CD4",
  "BUV395_CD11b",
  "BV650_Ly6G",
  "BV711_CD8a",
  "PECF594_CD80",
  "PECy7_NK1_1")

DrawNetworkPlot(
  dat = dat,
  timepoints = timepoints,
  timepoint.col = 'Group',
  cluster.col = 'TrackSOM_metacluster_lineage_tracking',
  marker.cols = markers,
  arrow.head.gap = 3,
  arrow.length = 3,
  img.height = 20,
  img.width = 20,
  min.node.size = 6,
  max.node.size = 12,
  load.temp.data = FALSE,
  standard.colours = "inferno",
  graph.layout = 'fr',
  calculation.type = 'mean',
  legend.position = 'right',
  file.format = "pdf",
  line.width = 0.8,
  no_merge = TRUE
)

setwd(work_dir)
dir.create("network_spectral")
setwd("network_spectral")

DrawNetworkPlot(
  dat = dat,
  timepoints = timepoints,
  timepoint.col = 'Group',
  cluster.col = 'TrackSOM_metacluster_lineage_tracking',
  marker.cols = c("FITC_Ly6C"),
  arrow.head.gap = 3,
  arrow.length = 3,
  img.height = 20,
  img.width = 20,
  min.node.size = 6,
  max.node.size = 12,
  load.temp.data = FALSE,
  standard.colours = "spectral",
  graph.layout = 'fr',
  calculation.type = 'mean',
  legend.position = 'none',
  file.format = "pdf",
  line.width = 0.8,
  no_merge = TRUE
)


