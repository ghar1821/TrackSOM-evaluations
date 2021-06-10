library(testthat)
library(data.table)

source("/Users/givanna/Documents/phd/code/FlowSOM-tracking/evaluation/calculate_tracking_accuracy.R")
source("/Users/givanna/Documents/phd/code/FlowSOM-tracking/evaluation/map_clusters.R")


test_row_found <- function(data, row) {
  return(nrow(merge(data, row)) > 0)
}

test_that("Edges are correctly inferred", {
  dataset <- data.table(
    SCA1=sample(seq(100,400), size=30),
    timepoint=c(
      rep(1,11),
      rep(2,14),
      rep(3,5)
    ),
    lineage=c(
      rep("A",3),
      rep("B",4),
      rep("C",4),
      rep("A",3),
      rep("A|1",3),
      rep("(B,C)",4),
      rep("D",4),
      rep("(A,A|1)",2),
      rep("(B,C),D)",3)
    ),
    association=c(
      rep("None",11),
      rep("A",6),
      rep("B & C",4),
      rep("A",4),
      rep("A & A|1",2),
      rep("(B,C) & D",3)
    )
  )
  
  # shuffle rows
  dataset <- dataset[sample(nrow(dataset)),]
  
  edges <- get_cluster_transitions(dat=dataset, 
                                   timepoints = c(1,2,3), 
                                   timepoint_col = "timepoint", 
                                   cluster_lineage_col = "lineage",
                                   cluster_assoc_col = "association")
  
  expected_edges <- data.frame(
    from_timepoint=c(rep(1,5),rep(2,4)),
    from=c("A","A","B","C","A","A","A|1","(B,C)","D"),
    to_timepoint=c(rep(2,5),rep(3,4)),
    to=c("A","A|1","(B,C)","(B,C)","D","(A,A|1)","(A,A|1)","(B,C),D)","(B,C),D)")
  )
  
  for (i in c(1:nrow(expected_edges))) {
    expect_true(test_row_found(edges, expected_edges[i,]))  
  }
})

test_that("Edges are correctly inferred when not all cells are clustered", {
  dataset <- data.table(
    SCA1=sample(seq(100,400), size=33),
    timepoint=c(
      rep(1,11),
      rep(2,14),
      rep(3,5),
      1,2,3
    ),
    lineage=c(
      rep("A",3),
      rep("B",4),
      rep("C",4),
      rep("A",3),
      rep("A|1",3),
      rep("(B,C)",4),
      rep("D",4),
      rep("(A,A|1)",2),
      rep("(B,C),D)",3),
      rep(NA,3)
    ),
    association=c(
      rep("None",11),
      rep("A",6),
      rep("B & C",4),
      rep("A",4),
      rep("A & A|1",2),
      rep("(B,C) & D",3),
      rep(NA,3)
    )
  )
  # shuffle rows
  dataset <- dataset[sample(nrow(dataset)),]
  edges <- get_cluster_transitions(dat=dataset,
                                   timepoints = c(1,2,3),
                                   timepoint_col = "timepoint",
                                   cluster_lineage_col = "lineage",
                                   cluster_assoc_col = "association")
  
  expected_edges <- data.frame(
    from_timepoint=c(rep(1,5),rep(2,4)),
    from=c("A","A","B","C","A","A","A|1","(B,C)","D"),
    to_timepoint=c(rep(2,5),rep(3,4)),
    to=c("A","A|1","(B,C)","(B,C)","D","(A,A|1)","(A,A|1)","(B,C),D)","(B,C),D)")
  )
  
  for (i in c(1:nrow(expected_edges))) {
    expect_true(test_row_found(edges, expected_edges[i,]))
  }
  
  expect_edges_not_exists <- data.frame(
    from_timepoint=c(1,2),
    to_timepoint=c(2,3),
    from=c("NA","NA"),
    to=c("NA","NA")
  )
  for (i in c(1:nrow(expect_edges_not_exists))) {
    expect_false(test_row_found(edges, expect_edges_not_exists[i,]))
  }
  
  expect_edges_not_exists <- data.frame(
    from_timepoint=c(1,2),
    to_timepoint=c(2,3),
    from=c(NA,NA),
    to=c(NA,NA)
  )
  for (i in c(1:nrow(expect_edges_not_exists))) {
    expect_false(test_row_found(edges, expect_edges_not_exists[i,]))
  }
  
})

test_that("Perfect transitions are correctly computed", {
  dataset <- data.table(
    SCA1=sample(seq(100,400), size=30),
    timepoint=c(
      rep(1,11),
      rep(2,14),
      rep(3,5)
    ),
    population=c(
      rep("T cell",3),
      rep("Monoblast",4),
      rep("Monocyte",4),
      rep("T cell",3),
      rep("T cell",3),
      rep("Monocyte",4),
      rep("T cell",4),
      rep("T cell",2),
      rep("Monocyte",2),
      rep("T cell",1)
    ),
    lineage=c(
      rep("A",3),
      rep("B",4),
      rep("C",4),
      rep("A",3),
      rep("A|1",3),
      rep("(B,C)",4),
      rep("D",4),
      rep("(A,A|1)",2),
      rep("(B,C)",2),
      rep("D",1)
    ),
    association=c(
      rep("None",11),
      rep("A",6),
      rep("B & C",4),
      rep("A",4),
      rep("A & A|1",2),
      rep("(B,C)",2),
      rep("D",1)
    )
  )
  
  rule = data.table(
    from=c("T cell", "Monoblast", "Monoblast", "Monocyte"),
    to=c("T cell", "Monoblast", "Monocyte", "Monocyte")
  )
  
  tracking_acc <- calculate_tracking_accuracy(dat = dataset,
                                              timepoints = c(1:3),
                                              timepoint_col = "timepoint",
                                              cluster_lineage_col = "lineage",
                                              cluster_assoc_col = "association",
                                              population_col = "population",
                                              rule = rule)
  expect_equal(tracking_acc, 1.0)
})

test_that("Perfect transitions are correctly computed even with some cells not assigned", {
  dataset <- data.table(
    SCA1=sample(seq(100,400), size=33),
    timepoint=c(
      rep(1,11),
      rep(2,14),
      rep(3,5),
      1,2,3
    ),
    population=c(
      rep("T cell",3),
      rep("Monoblast",4),
      rep("Monocyte",4),
      rep("T cell",3),
      rep("T cell",3),
      rep("Monocyte",4),
      rep("T cell",4),
      rep("T cell",2),
      rep("Monocyte",2),
      rep("T cell",1),
      rep("B cell",3)
    ),
    lineage=c(
      rep("A",3),
      rep("B",4),
      rep("C",4),
      rep("A",3),
      rep("A|1",3),
      rep("(B,C)",4),
      rep("D",4),
      rep("(A,A|1)",2),
      rep("(B,C)",2),
      rep("D",1),
      rep(NA,3)
    ),
    association=c(
      rep("None",11),
      rep("A",6),
      rep("B & C",4),
      rep("A",4),
      rep("A & A|1",2),
      rep("(B,C)",2),
      rep("D",1),
      rep(NA,3)
    )
  )
  
  rule = data.table(
    from=c("T cell", "Monoblast", "Monoblast", "Monocyte"),
    to=c("T cell", "Monoblast", "Monocyte", "Monocyte")
  )
  
  tracking_acc <- calculate_tracking_accuracy(dat = dataset,
                                              timepoints = c(1:3),
                                              timepoint_col = "timepoint",
                                              cluster_lineage_col = "lineage",
                                              cluster_assoc_col = "association",
                                              population_col = "population",
                                              rule = rule)
  expect_equal(tracking_acc, 1.0)
})

test_that("Some illegal transitions in tracking by lineage", {
  dataset <- data.table(
    SCA1=sample(seq(100,400), size=10),
    timepoint=c(
      rep(1,5),
      rep(2,5)
    ),
    population=c(rep('T cell',3),
                 rep('B cell', 2),
                 rep('T cell', 5)
    ),
    lineage=c(rep('A',3),
              rep('B', 2),
              rep('(A,B)', 5)
    ),
    association=c(rep('None',5),
                  rep('A & B',5)
    )
  )
  
  rule = data.table(
    from=c("T cell", "B cell"),
    to=c("T cell", "B cell")
  )
  
  tracking_acc <- calculate_tracking_accuracy(dat = dataset,
                                              timepoints = c(1:2),
                                              timepoint_col = "timepoint",
                                              cluster_lineage_col = "lineage",
                                              cluster_assoc_col = "association",
                                              population_col = "population",
                                              rule = rule)
  expect_equal(tracking_acc, 0.5)
})

test_that("Some illegal transitions in tracking by association", {
  dataset <- data.table(
    SCA1=sample(seq(100,400), size=15),
    timepoint=c(
      rep(1,5),
      rep(2,10)
    ),
    population=c(rep('T cell',3),
                 rep('T cell', 2),
                 rep('T cell', 5),
                 rep('T cell', 2),
                 rep('Eosinophil', 3)
    ),
    lineage=c(rep('A',3),
              rep('B', 2),
              rep('(A,B)', 5),
              rep('B|1',2),
              rep('C', 3)
    ),
    association=c(rep('None',5),
                  rep('A & B',5),
                  rep('B', 2),
                  rep('A',3)
    )
  )
  
  rule = data.table(
    from=c("T cell", "B cell", 'Eosinophil'),
    to=c("T cell", "B cell", 'Eosinophil')
  )
  
  tracking_acc <- calculate_tracking_accuracy(dat = dataset,
                                              timepoints = c(1:2),
                                              timepoint_col = "timepoint",
                                              cluster_lineage_col = "lineage",
                                              cluster_assoc_col = "association",
                                              population_col = "population",
                                              rule = rule)
  expect_equal(tracking_acc, 0.75)
})

test_that("Some illegal transitions in tracking by lineage and association", {
  dataset <- data.table(
    SCA1=sample(seq(100,400), size=15),
    timepoint=c(
      rep(1,5),
      rep(2,10)
    ),
    population=c(rep('T cell',3),
                 rep('B cell', 2),
                 rep('T cell', 5),
                 rep('B cell', 2),
                 rep('Eosinophil', 3)
    ),
    lineage=c(rep('A',3),
              rep('B', 2),
              rep('(A,B)', 5),
              rep('B|1',2),
              rep('C', 3)
    ),
    association=c(rep('None',5),
                  rep('A & B',5),
                  rep('B', 2),
                  rep('A',3)
    )
  )
  
  rule = data.table(
    from=c("T cell", "B cell", 'Eosinophil'),
    to=c("T cell", "B cell", 'Eosinophil')
  )
  
  tracking_acc <- calculate_tracking_accuracy(dat = dataset,
                                              timepoints = c(1:2),
                                              timepoint_col = "timepoint",
                                              cluster_lineage_col = "lineage",
                                              cluster_assoc_col = "association",
                                              population_col = "population",
                                              rule = rule)
  expect_equal(tracking_acc, 0.5)
})