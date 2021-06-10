library(data.table)
library(Spectre)
library(stringr)
library(scales)
library(TrackSOM)

work_dir <- "~/Dropbox (Sydney Uni)/tracksom/wnv_cns/clustered_10x10_pv/"
setwd(work_dir)

dat <- fread("res_10x10_pv.csv")

markers <- names(dat)[c(1:8,10:20)]
timepoints <- c("Mock","WNV-01","WNV-02","WNV-03","WNV-04","WNV-05")

dat_sub <- do.subsample(dat,
                        targets = rep(50000, 6),
                        divide.by = 'Group')
dat_sub <- run.fitsne(dat_sub, use.cols = markers)

setwd(work_dir)
setwd("fitsne_map/")
dir.create("density_plot")
setwd("density_plot")

n_cell_per_grp <- dat[, .(cnt_total=.N), by='Group']
labels <- sapply(c(1: nrow(n_cell_per_grp)), function(i) {
    row <- n_cell_per_grp[i,]
    grp <- row$Group
    paste(grp, "\n50,000 out of ", 
          prettyNum(row$cnt_total, big.mark = ',', scientific=FALSE),
          "cells"
    )
})
names(labels) <- n_cell_per_grp$Group

# Spectre's make.colour.plot with some modifications
p <- ggplot(data = dat_sub, aes(x = .data[["FItSNE_X"]], y = .data[["FItSNE_Y"]])) + 
    ggpointdensity::geom_pointdensity(size = 1) +
    ggplot2::scale_colour_gradientn(
        colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(50))) +
    facet_wrap( ~ Group, labeller = labeller(Group = labels)) +
    theme_bw() + 
    labs(color='Density')
ggsave("density_facet.png", p,
       width = 10,
       height = 7,
       dpi = 800,
       limitsize = FALSE)
