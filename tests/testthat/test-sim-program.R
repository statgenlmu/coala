context('Simulation Program Creation')

test_that("test.create_simprog", {
    expect_true(exists("simprograms"))
    simprog_nr <- length(simprograms)
    #expect_error(create_simprog("name", "feature", sin, cos))
    #expect_error(create_simprog("name"))

    create_simprog("test1", "feature", "sum.stat", sin)
    expect_equal(length(simprograms), simprog_nr + 1)

    create_simprog("test2", "feature", "sum.stat", sin, cos)
    expect_equal(length(simprograms), simprog_nr + 2)

    create_simprog("test2", "feature", "sum.stat", sin, cos, tan)
    expect_equal(length(simprograms), simprog_nr + 2)
})
