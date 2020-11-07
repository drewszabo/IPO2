test_that("%nin% works", {
  expect_equal((c(1, 2, 3) %nin% c(3, 4, 5)), c(TRUE, TRUE, FALSE))
})
