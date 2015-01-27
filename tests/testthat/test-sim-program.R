context('Simulation Program Creation')

test_that("test.createSimProgram", {
    expect_true(exists("sim_programs"))
    sim_prog_nr = length(sim_programs)
    #expect_error(createSimProgram("name", "feature", sin, cos))
    #expect_error(createSimProgram("name"))

    createSimProgram("test1", "feature", "sum.stat", sin)
    expect_equal(length(sim_programs), sim_prog_nr + 1)

    createSimProgram("test2", "feature", "sum.stat", sin, cos)
    expect_equal(length(sim_programs), sim_prog_nr + 2)

    createSimProgram("test2", "feature", "sum.stat", sin, cos, tan)
    expect_equal(length(sim_programs), sim_prog_nr + 2)
})

