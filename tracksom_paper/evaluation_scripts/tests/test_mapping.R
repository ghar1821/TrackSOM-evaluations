library(testthat)
library(data.table)

source("/Users/givanna/Documents/phd/code/FlowSOM-tracking/evaluation/map_clusters.R")

test_row_found <- function(data, row) {
  return(nrow(merge(data, row)) > 0)
}

test_that("Perfect clusters 1 time point", {
  # 2 populations (T cell and Monocyte) split perfectly into 2 clusters (A and B)
  dataset <- data.table(
    SCA1=sample(400, 10),
    population=c(rep("T cell", 6),
                 rep("Monocyte", 4)),
    cluster=c(rep("A", 6),
              rep("B", 4)),
    timepoint=c(rep(1,10))
  )
  # shuffle data
  dataset <- dataset[sample(nrow(dataset)),]
  
  mapping <- data.frame(map_clusters(dataset, "timepoint", 'cluster', 'population'))
  
  expect_equal(mapping[mapping$cluster == 'A', c(2)], 'T cell')
  expect_equal(mapping[mapping$cluster == 'B', c(2)], 'Monocyte')
})

test_that("Perfect clusters 2 time points", {
  # 2 populations (T cell and Monocyte) split perfectly into 2 clusters (A and B)
  dataset <- data.table(
    SCA1=sample(400, 20),
    population=c(rep("T cell", 6),
                 rep("Monocyte", 4),
                 rep("T cell", 6),
                 rep("Monocyte", 4)),
    cluster=c(rep("A", 6),
              rep("B", 4),
              rep("A", 6),
              rep("B", 4)),
    timepoint=c(rep(1,10),
                rep(2,10))
  )
  # shuffle data
  dataset <- dataset[sample(nrow(dataset)),]
  
  mapping <- data.frame(map_clusters(dataset, "timepoint", 'cluster', 'population'))
  expect_equal(mapping[mapping$cluster == 'A' & mapping$timepoint == 1, c(2)], 'T cell')
  expect_equal(mapping[mapping$cluster == 'B' & mapping$timepoint == 1, c(2)], 'Monocyte')
  expect_equal(mapping[mapping$cluster == 'A' & mapping$timepoint == 2, c(2)], 'T cell')
  expect_equal(mapping[mapping$cluster == 'B' & mapping$timepoint == 2, c(2)], 'Monocyte')
  
})

test_that("Populations split between 2 clusters", {
  # 2 populations (T cell and Monocyte) split into 4 clusters (A-D)
  dataset <- data.table(
    SCA1=sample(400, 20),
    population=c(rep("T cell", 8),
                 rep("Monocyte", 4),
                 rep("T cell", 4),
                 rep("Monocyte", 4)),
    cluster=c(rep("A", 4),
              rep("B", 4),
              rep("C", 4),
              rep("A|1", 1),
              rep("(A,B)", 3),
              rep("C", 4)),
    timepoint=c(rep(1,12),
                rep(2,8))
  )
  # shuffle data
  dataset <- dataset[sample(nrow(dataset)),]

  mapping <- data.frame(map_clusters(dataset, "timepoint", 'cluster', 'population'))
  
  expect_equal(mapping[mapping$cluster == 'A' & mapping$timepoint == 1, c(2)], 'T cell')
  expect_equal(mapping[mapping$cluster == 'B' & mapping$timepoint == 1, c(2)], 'T cell')
  expect_equal(mapping[mapping$cluster == 'C' & mapping$timepoint == 1, c(2)], 'Monocyte')
  expect_equal(mapping[mapping$cluster == 'A|1' & mapping$timepoint == 2, c(2)], 'T cell')
  expect_equal(mapping[mapping$cluster == '(A,B)' & mapping$timepoint == 2, c(2)], 'T cell')
  expect_equal(mapping[mapping$cluster == 'C' & mapping$timepoint == 2, c(2)], 'Monocyte')
})

test_that("Clusters containing multiple populations", {
  # 2 populations (T cell and Monocyte) split into 3 clusters with 1 cluster containing
  # mix of T cell and Monocyte.
  dataset <- data.table(
    SCA1=sample(400, 7),
    population=c(rep("T cell", 3),
                 rep("Monocyte", 4)),
    cluster=c('A','B','B','A','C','C','C'),
    timepoint=c(rep(1,7))
  )
  # shuffle data
  dataset <- dataset[sample(nrow(dataset)),]

  mapping <- data.frame(map_clusters(dataset, "timepoint", 'cluster', 'population'))
  
  expect_true(mapping[mapping$cluster == 'A' & mapping$timepoint == 1, c(2)] == 'T cell' ||
                mapping[mapping$cluster == 'A' & mapping$timepoint == 1, c(2)] == 'Monocyte')
  expect_equal(mapping[mapping$cluster == 'B' & mapping$timepoint == 1, c(2)], 'T cell')
  expect_equal(mapping[mapping$cluster == 'C' & mapping$timepoint == 1, c(2)], 'Monocyte')
  
  dataset <- data.table(
    SCA1=sample(400, 7),
    population=c(rep("Monocyte", 4),
                 rep("T cell", 3)),
    cluster=c('C','C','C','A','B','B','A'),
    timepoint=c(rep(1,7))
  )
})

# test_that("Some cells are not assigned any cluster", {
#   # 2 populations (T cell and Monocyte),Monocyte not assigned any cluster
#   dataset <- data.table(
#     SCA1=sample(400, 3),
#     population=c(rep("T cell", 2),
#                  rep("Monocyte", 1)),
#     cluster=c('A','A',NA),
#     timepoint=rep(1,3)
#   )
#   # shuffle data
#   dataset <- dataset[sample(nrow(dataset)),]
# 
#   mapping <- data.frame(map_clusters(dataset, "timepoint", 'cluster', 'population'))
#   expect_equal(mapping[mapping$cluster == 'A' & mapping$timepoint == 1, c(2)], 'T cell')
# })