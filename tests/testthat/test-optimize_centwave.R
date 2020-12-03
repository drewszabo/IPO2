test_that("check_centwave_params works", {
  expect_error(check_centwave_params(list(foo = "bar")), "These parameters are not used by centWave: ")
  expect_error(check_centwave_params(list(ppm = c(1, 2, 3))), "Too many parameters specified for: ")
  expect_error(check_centwave_params(list(ppm = 1)), "No parameters specified for optimization")
  expect_error(check_centwave_params(list(ppm = c(-1, 1))), "must be greater than 0")
  expect_error(check_centwave_params(list(snthresh = c(-1, 1))), "must not be negative")
})
