ugene <- inferNetwork(Repressilator, mtry=3L)
test_that("cutoff too high", {
  expect_error(tuneThreshold(Repressilator, ugene, cutoffs = c(0.01, 0.99)))
})

# runs without errors
result <- tuneThreshold(Repressilator, ugene)
test_that("return values", {
  expect_equal(class(result), "ugene.analysis")
  expect_equal(length(result$step.errors), 31)
  expect_equal(length(result$col.errors), 5)
  expect_equal(length(result$step.masks), 31)
  expect_equal(length(result$col.masks), 5)
  expect_equal(dim(result$step.masks[[1]])[1], 6)
  expect_equal(dim(result$step.masks[[1]])[2], 6)
  expect_equal(dim(result$col.masks[[1]])[1], 6)
  expect_equal(dim(result$col.masks[[1]])[2], 6)
})
