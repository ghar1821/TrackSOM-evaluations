library(TrackSOM)
library(Spectre)
library(tools)

data_dir <- "~/Documents/phd/wnv_cns/data/"
setwd(data_dir)
dat_per_tp <- read.files()
dat <- rbindlist(dat_per_tp)

meta_dat <- fread("~/Documents/phd/wnv_cns/metadata/sample.details.csv")
setnames(meta_dat, 'Filename', 'FileName')
dat <- do.add.cols(dat, base.col = 'FileName', add.dat = meta_dat, add.by = 'FileName', rmv.ext = TRUE)

groups <- c('Mock', 'WNV-01', 'WNV-02', 'WNV-03', 'WNV-04', 'WNV-05')

# Make sure this is sorted based on group above
dat_per_tp <- split(dat, dat$Group)

markers <- names(dat)[c(1:8, 10:12, 14:20)]

message("Running TrackSOM")
tracksom.result <- TrackSOM(inputFiles = dat_per_tp,
                            colsToUse = markers,
                            tracking = TRUE,
                            noMerge = TRUE,  # TODO change me if you want to allow merging
                            seed = 42,
                            xdim = 10,
                            ydim = 10,
                            nClus = c(25, 27, 30, 30, 30, 30),
                            scale = TRUE,
                            dataFileType = "data.frame"  # TODO change me according to file type you have
)
dat_clustered <- ConcatenateClusteringDetails(tracksom.result = tracksom.result,
                                             dat = dat,
                                             timepoint.col = 'Group',
                                             timepoints = groups)

save_dir <- "~/Documents/phd/wnv_cns/clustered_10x10_PV"
dir.create(save_dir)
setwd(save_dir)
fwrite(dat_clustered, "res_10x10_PV.csv")
write.files(dat_clustered, "res_10x10_PV", divide.by = 'Group')

