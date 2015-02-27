context('sumstat File')

test_that('File statistic works', {
  folder <- tempfile('sumstat_file_test')
  stat <- sumstat_file(folder, 3)

  files <- c(tempfile("test1"), tempfile("test2"))
  cat('test1', file = files[1])
  cat('test2', file = files[2])

  expect_equal(stat$get_name(), 'file')
  expect_equal(stat$get_group(), 3)
  files_copy <- stat$calculate(NULL, files, NULL)

  expect_true(file.exists(folder))
  expect_true(all(file.exists(files_copy)))
  expect_equal(scan(files_copy[1], what='character', quiet = TRUE), 'test1')
  expect_equal(scan(files_copy[2], what='character', quiet = TRUE), 'test2')

  unlink(c(files, folder), recursive = TRUE)
})
