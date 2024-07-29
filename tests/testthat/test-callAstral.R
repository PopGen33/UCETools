test_that("Test that ASTRAL can be downloaded and called", {
    # Implicitly tests getAstral() as well
    # Only tests that it can be called by asking for the help message
    temp_dir <- tempdir()
    getAstral(install.dir = temp_dir)
    expect_true(file.exists(normalizePath(file.path(temp_dir, "ASTRAL-5.7.1", "Astral", "astral.5.7.1.jar"))))
    expect_equal(normalizePath(getOption("astralPath")), normalizePath(file.path(temp_dir, "ASTRAL-5.7.1", "Astral", "astral.5.7.1.jar")))

    # Check for that callAstral returns 1 when asking for the help message
    # Why ASTRAL returns exit code 1 when successfully asking for the help message is beyond me
    # This is, therefore, a terrible test
    expect_equal(callAstral(help = TRUE), 1)
})
