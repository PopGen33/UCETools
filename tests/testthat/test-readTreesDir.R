test_that("Tests that readTreesDir() reads the example trees correctly", {
    # R doesn't like tests running on more than 2 cores?
    ex_trees <- readTreesDir(system.file("extdata", "UCE_aligns_genetrees", package = "UCETools"), threads = 2)
    expect_equal(length(ex_trees), 800)
})
