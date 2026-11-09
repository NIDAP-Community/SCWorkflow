local({
	tests <- if (file.exists(file.path("..", "setup-fixtures.R"))) {
		".."
	} else {
		file.path(testthat::test_path(), "..")
	}
	source(file.path(tests, "setup-fixtures.R"), local = TRUE)
	scworkflow_setup_fixtures(root = file.path(tests, ".."),
														fixtures = file.path(tests, "testthat", "fixtures"))
})