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

setwd(work_dir)
mapping <- fread("meta_pop_mapping.csv")

meta_incl <- sapply(unique(dat$TrackSOM_metacluster_lineage_tracking), function(d) {
    return(str_extract_all(d, "[a-zA-Z]+")[[1]] %in% mapping$Metacluster)
})
meta_incl <- data.table(metacluster=names(meta_incl), include=meta_incl)
meta_incl <- meta_incl[meta_incl$include == TRUE]
meta_incl$origin <- sapply(meta_incl$metacluster, function(x) str_extract_all(x, "[a-zA-Z]+")[[1]])
meta_incl <- merge(meta_incl, mapping, by.x = 'origin', by.y = 'Metacluster')

dat_some_pop <- merge(dat, meta_incl, 
                      by.x = 'TrackSOM_metacluster_lineage_tracking',
                      by.y = 'metacluster',
                      all = TRUE)

# I|2 captured some neutrophils
dat_some_pop[dat_some_pop$TrackSOM_metacluster_lineage_tracking == 'I|2', Population := NA]

# this is the scatter plot per sample ----
cnt_per_pop <- dat_some_pop[, .(cnt_pop_grp = .N), by=c("Population", "Group", "Sample")]
cnt_per_grp <- dat_some_pop[, .(cnt_grp = .N), by=c("Sample", "Group")]
cnt_all <- merge(cnt_per_grp, cnt_per_pop)
cnt_all <- cnt_all[!is.na(cnt_all$Population),]
cnt_all[, proportion:=cnt_pop_grp/cnt_grp]

pops <- unique(cnt_all$Population)
for (p in pops[c(1,2,3,6)]) {
    cnt_all_p <- cnt_all[cnt_all$Population == p,]
    make.autograph(
        dat = cnt_all_p,
        x.axis = "Group",
        y.axis = "proportion",
        y.axis.label = "Proportion of cells",
        x.axis.label = "Disease severity",
        colour.by = "Group",
        colours = c("yellow1", "gold1", "orange2", "darkorange2", "red1", "red4"),
        title = p,
        filename = paste0("proportion_", p, ".pdf"),
        violin = FALSE,
        dot.size = 4,
        max.y = 1.1,
        max.y.value = NULL
    )
}

for (p in pops[c(4,5)]) {
    cnt_all_p <- cnt_all[cnt_all$Population == p,]
    make.autograph(
        dat = cnt_all_p,
        x.axis = "Group",
        y.axis = "proportion",
        y.axis.label = "Proportion of cells",
        x.axis.label = "Disease severity",
        colour.by = "Group",
        colours = c("yellow1", "gold1", "orange2", "darkorange2", "red1", "red4"),
        title = p,
        filename = paste0("proportion_", p, ".pdf"),
        violin = FALSE,
        dot.size = 4,
        max.y = 1.0,
        max.y.value = 1.0
    )
}

# draw count divided by number of mice
cnt_sample <- unique(cnt_all[, c("Sample", "Group")])
cnt_sample <- cnt_sample[, .(cnt_samp = .N), by='Group']
cnt_all <- merge(cnt_all, cnt_sample, by='Group')
cnt_all[, abs_cnt_div := cnt_pop_grp / cnt_samp]

for (p in pops) {
    cnt_all_p <- cnt_all[cnt_all$Population == p,]
    make.autograph(
        dat = cnt_all_p,
        x.axis = "Group",
        y.axis = "abs_cnt_div",
        y.axis.label = "Count of cells divided by number of mice",
        x.axis.label = "Disease severity",
        colour.by = "Group",
        colours = c("yellow1", "gold1", "orange2", "darkorange2", "red1", "red4"),
        title = p,
        filename = paste0("count_", p, ".pdf"),
        violin = FALSE,
        dot.size = 4,
        max.y = 1.1,
        max.y.value = NULL
    )
}

