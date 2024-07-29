# Setup for testing for the UCETools package; uses usethis()

# Load testing packages
library(testthat)
library(usethis)
library(devtools)
library(here)

#use_testthat()
#UCETools_simTree <- read.tree()
#use_data(UCETools_simTree, version = 3)
#UCETools_examplePISMatrix <- getPisCounts(here("tests", "UCE_aligns"))
#UCETools_examplePISMatrix[,1] <- NA
#use_data(UCETools_examplePISMatrix, version = 3)

trees_dir_ex <- system.file("extdata", "UCE_aligns_genetrees", package = "UCETools")
UCETools_genetrees <- readTreesDir(trees_dir_ex)
