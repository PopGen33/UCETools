test_that("Test that getPisCounts returns correctly", {
    # R doesn't like tests running on more than 2 cores?
    # Load alignments
    alignmentDir <- system.file("extdata", "UCE_aligns", package = "UCETools")
    # Check equality of useAmbig True/False; should be equal on this dataset
    ex_PISCounts <- getPisCounts(alignmentDir, threads = 2)
    ex_PISCounts_amb <- getPisCounts(alignmentDir, useAmbiguity = TRUE, threads = 2)
    expect_equal(ex_PISCounts,ex_PISCounts_amb)
    # Check summary stats of example data are equal to expectations
    expect_equal(min(unlist(ex_PISCounts[,3])), 0)
    expect_equal(max(unlist(ex_PISCounts[,3])), 108)
    expect_equal(median(unlist(ex_PISCounts[,3])), 8)
    expect_equal(mean(unlist(ex_PISCounts[,3])), 10.665)
    rm(alignmentDir, ex_PISCounts, ex_PISCounts_amb)
})
