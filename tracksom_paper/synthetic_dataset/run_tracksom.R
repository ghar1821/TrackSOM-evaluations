# The script was run on a HPC platform where a PBS script equipped with
# the hyperparameter values were sourced.

# Import libraries
library(TrackSOM)
library(Spectre)

# CHANGE THIS

setwd("/project/chronoclust/jon/FlowSOM-tracking/synthetic_dataset/dataset/data_coarse_noNoise/fcs")
PrimaryDirectory <- getwd()


data.files <- list.files(PrimaryDirectory, ".fcs")
# convert to absolute path
data.files.fullpath <- sapply(data.files, function(fname) {
    full.fname <- paste(PrimaryDirectory, fname, sep="/")
})

# get command line args. These can be read in from LHC samples
args <- commandArgs(trailingOnly = TRUE)
meta.max <- round(as.numeric(args[1]))
meta.per.timepoint <- as.list(c(round(as.numeric(args[2])), round(as.numeric(args[3])), round(as.numeric(args[4])), round(as.numeric(args[5])), round(as.numeric(args[6]))))
grid.size <- round(as.numeric(args[7]))
option <- as.integer(args[8])
outdir <- args[9]

# column to be used by TrackSOM to do clustering
ClusteringCols <- c("x", "y", "z")

## AA
if (option == 1) {
    
    tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                                colsToUse = ClusteringCols,
                                tracking = TRUE,
                                seed = 42,
                                xdim = grid.size,
                                ydim = grid.size,
                                maxMeta = meta.max
    )
} else if (option == 2) {
    ## PI
    ## Option 2: FlowSOM creates same number of meta clusters per time point.
    tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                                colsToUse = ClusteringCols,
                                tracking = TRUE,
                                seed = 42,
                                xdim = grid.size,
                                ydim = grid.size,
                                nClus = meta.max
    )
} else if (option == 3) {
    ## PV
    ## Option 3: FlowSOM creates different number of meta clusters per time point.
    print(meta.per.timepoint)
    tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                                colsToUse = ClusteringCols,
                                tracking = TRUE,
                                seed = 42,
                                xdim = grid.size,
                                ydim = grid.size,
                                nClus = meta.per.timepoint
    )
}

cell.dat <- Spectre::read.files(file.loc = PrimaryDirectory,
                                file.type = ".fcs")
cell.dat <- Spectre::do.merge.files(cell.dat)
head(cell.dat)

cell.dat <- consolidate.tracksom.result(tracksom.result = tracksom.result,
                                        dat = cell.dat,
                                        divide.by = "FileName")
setwd(outdir)
Spectre::write.files(cell.dat, "Result", divide.by = "FileName")