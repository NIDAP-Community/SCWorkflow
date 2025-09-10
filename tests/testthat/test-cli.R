
test_that("scworkflow helps", {
  expect_snapshot(cli_exec("--help"))
  expect_snapshot(system(paste(system.file("exec", "scworkflow", package = "SCWorkflow"), "--help")))
  expect_true(inherits(cli_exec(c("filterQC", "--help")), "help_files_with_topic"))
  expect_warning(cli_exec("not_a_function"), "not a known function")
})