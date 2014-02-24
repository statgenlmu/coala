context('Simulation programs')

test_that('list of simulation programs exists', {
  expect_is(sim.progs, 'list')
  expect_is(sim.progs$ms, 'environment')
})

test_that('all simulation programs are complete', {
  for (sim.prog in sim.progs) {
    expect_is(sim.prog$simulate, 'function')
    expect_is(sim.prog$finalize, 'function')
    expect_is(sim.prog$features, 'character')
    expect_is(sim.prog$sum.stats, 'character')
  }
})
