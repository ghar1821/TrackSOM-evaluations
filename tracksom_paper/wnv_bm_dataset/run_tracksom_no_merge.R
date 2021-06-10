# Import libraries
library(TrackSOM)
library(Spectre)
library(data.table)

# CHANGE THIS

setwd("/project/chronoclust/tracksom_eval/tracksom_no_merge/wnv_dataset")
PrimaryDirectory <- getwd()


data.files <- list.files(PrimaryDirectory, ".fcs")
# convert to absolute path
data.files.fullpath <- sapply(data.files, function(fname) {
    full.fname <- paste(PrimaryDirectory, fname, sep="/")
})

# get command line args
# This came from the PBS script in the HPC platform. Read them in from LHC samples
args <- commandArgs(trailingOnly = TRUE)
meta.per.timepoint <- as.list(c(round(as.numeric(args[2])), round(as.numeric(args[3])), round(as.numeric(args[4])), round(as.numeric(args[5])), round(as.numeric(args[6])), round(as.numeric(args[7])), round(as.numeric(args[8])), round(as.numeric(args[9]))))
grid.size <- round(as.numeric(args[10]))
option <- as.integer(args[11])
outdir <- args[12]
meta.max <- round(as.numeric(args[1]))

ClusteringCols <- c("FSC.A.", "BV711.CD11c.", "BV785.B220.", "PE.CD115.", "PECy5.CD3.CD19.NK11.", "PECy7.CD16.32.", "SSC.A.", "FITC.Ly6C.", "PerCPCy5.5.CD45.", "AF700.CD48.", "APCCy7.Ly6G.", "BV421.CD117.", "BV510.SCA.1.", "BV605.CD11b.")
print(option)
if (option == 1) {
    ## Option 1: FlowSOM infer the optimal number of meta clusters. 
    ## AA
    tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                                colsToUse = ClusteringCols,
                                tracking = TRUE,
                                seed = 42,
                                xdim = grid.size,
                                ydim = grid.size,
                                #rlen = r.len,
                                #alpha = alpha.range,
                                maxMeta = meta.max,
                                noMerge = TRUE,
                                dataFileType = ".fcs"
    )
} else if (option == 2) {
    ## Option 2: FlowSOM creates same number of meta clusters per time point.
    ## PI
    tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                                colsToUse = ClusteringCols,
                                tracking = TRUE,
                                seed = 42,
                                xdim = grid.size,
                                ydim = grid.size,
                                #rlen = r.len,
                                #alpha = alpha.range,
                                nClus = meta.max,
                                noMerge = TRUE,
                                dataFileType = ".fcs"
    )
} else if (option == 3) {
    ## Option 3: FlowSOM creates different number of meta clusters per time point.
    ## PV
    tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                                colsToUse = ClusteringCols,
                                tracking = TRUE,
                                seed = 42,
                                xdim = grid.size,
                                ydim = grid.size,
                                nClus = meta.per.timepoint,
                                noMerge = TRUE,
                                dataFileType = ".fcs"
    )
}

cell.dat <- Spectre::read.files(file.loc = PrimaryDirectory,
                                file.type = ".fcs")
cell.dat <- Spectre::do.merge.files(cell.dat)
head(cell.dat)

cell.dat <- ConcatenateClusteringDetails(tracksom.result = tracksom.result,
                                         dat = cell.dat,
                                         timepoint.col = "FileName",
                                         timepoints = sapply(c(0:7), function(d) paste0("WNV_D", d)))
setwd(outdir)
Spectre::write.files(cell.dat, "Result", divide.by = "FileName")