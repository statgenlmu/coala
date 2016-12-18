context("Simulators")

progs <- list(scrm = get_simulator("scrm"))
if (has_ms()) progs[["ms"]] <- get_simulator("ms")
if (has_msms()) progs[["msms"]] <- get_simulator("msms")

for (simulator in progs) {
  name <- simulator$get_name()

  test_that(paste(name, "can create simulation tasks"), {
    model <- coal_model(10, 1, 100) + feat_mutation(5) + sumstat_seg_sites()
    task <- simulator$create_task(model, NULL, 1)
    expect_equal(task$locus_number, 1)

    model <- coal_model(10, 2, 100) + feat_mutation(5) + sumstat_seg_sites()
    task <- simulator$create_task(model, NULL, 2)
    expect_equal(task$locus_number, 2)
  })


  test_that(paste(name, "can simulate segsites"), {
    model <- coal_model(10, 2, 100) + feat_mutation(5) + sumstat_seg_sites()
    task <- simulator$create_task(model, NULL, 2)
    result <- simulator$simulate(model, task)

    expect_that(result$seg_sites, is_a("list"))
    expect_equal(length(result$seg_sites), 2)
    expect_true(is_segsites(result$seg_sites[[1]]))
    expect_true(is_segsites(result$seg_sites[[2]]))
    expect_equal(nrow(result$seg_sites[[1]]), 10)
    expect_equal(nrow(result$seg_sites[[2]]), 10)
    expect_identical(result$simulator, simulator)
  })


  test_that(paste(name, "can simulate trees"), {
    model <- coal_model(10, 2, 100) + sumstat_trees()
    task <- simulator$create_task(model, NULL, 2)
    result <- simulator$simulate(model, task)

    expect_that(result$trees, is_a("list"))
    expect_equal(length(result$trees), 2)
    expect_that(result$trees[[1]], is_a("character"))
    expect_that(result$trees[[2]], is_a("character"))
    expect_equal(length(result$trees[[1]]), 1)
    expect_equal(length(result$trees[[2]]), 1)


    model <- coal_model(10, 2, 100) + feat_recombination(5) + sumstat_trees()
    task <- simulator$create_task(model, NULL, 2)
    result <- simulator$simulate(model, task)

    expect_equal(length(result$trees), 2)
    expect_that(result$trees[[1]], is_a("character"))
    expect_that(result$trees[[2]], is_a("character"))
    expect_gt(length(result$trees[[1]]), 1)
    expect_gt(length(result$trees[[2]]), 1)
  })


  test_that(paste(name, "can simulate files"), {
    model <- coal_model(10, 2, 100) +
      feat_mutation(5) +
      sumstat_file(".")
    task <- simulator$create_task(model, NULL, 2)
    result <- simulator$simulate(model, task)

    expect_true(!is.null(result$file))
    expect_true(file.exists(result$file))

    unlink(result$file)
  })


  test_that(paste(name, "simulations are reproducible"), {
    model <- coal_model(10, 3, 100) +
      feat_mutation(5) +
      sumstat_trees() +
      sumstat_seg_sites()
    task <- simulator$create_task(model, NULL, 3)

    set.seed(17); result1 <- simulator$simulate(model, task)
    set.seed(17); result2 <- simulator$simulate(model, task)
    expect_equal(result1, result2)
  })


  test_that(paste(name, "returns the simulation cmds"), {
    model <- coal_model(3, 2, 100) + sumstat_trees()
    task <- simulator$create_task(model, NULL, 2)
    result <- simulator$simulate(model, task)
    expect_true(grepl("3 2 -T", result$cmd))
  })


  test_that(paste(name, "can print simulation cmds"), {
    model <- coal_model(3, 2, 100) +
      feat_mutation(par_named("theta")) +
      sumstat_sfs()
    cmd <- simulator$get_cmd(model)
    expect_true(grepl(name, cmd))
    expect_true(grepl("3 2", cmd))
    expect_true(grepl("-t theta", cmd))

    model <- coal_model(2:3, 7, 100) +
      feat_mutation(par_named("theta")) +
      sumstat_sfs()
    cmd <- simulator$get_cmd(model)
    expect_true(grepl(name, cmd))
    expect_true(grepl("5 7", cmd))
    expect_true(grepl("-t theta", cmd))
  })


  test_that("unknown features generate an error", {
    feat_null <- R6Class("NULL", inherit = feature_class)$new()
    expect_error(simulator$create_cmd_template(model + feat_null))
  })
}
