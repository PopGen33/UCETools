test_that("Tests that readTreesDir() reads the example trees correctly", {
  ex_trees <- readTreesDir(system.file("extdata", "UCE_aligns_genetrees", package = "UCETools"))
  expect_equal(length(ex_trees), 800)
})
