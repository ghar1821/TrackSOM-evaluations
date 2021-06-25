# Remove the noise data points from the synthetic dataset as we don't need it.

library(Spectre)
library(Biobase)

setwd("~/Documents/phd/code/FlowSOM-tracking/synthetic_dataset/dataset")

data_list <- read.files(paste0(getwd(), '/gating_fine'))
data <- do.merge.files(data_list)

data <- data[data$PopName != 'Noise']

dir.create("data_fine_noNoise")
setwd("data_fine_noNoise")
write.files(dat = data,
            file.prefix = "synthetic",
            divide.by = 'FileName')

timepoints <- unique(data$FileName)
for (t in timepoints) {
  data_markers_only <- data[data$FileName == t, c(1:3)]
  write.files(dat = data_markers_only,
              file.prefix = t,
              write.fcs = TRUE)
}

test_files <- read.files(
  '~/Documents/phd/code/FlowSOM-tracking/synthetic_dataset/dataset/data_fine_noNoise/fcs', 
  file.type = ".fcs")
test_data <- do.merge.files(test_files)
head(test_data)
head(data)
