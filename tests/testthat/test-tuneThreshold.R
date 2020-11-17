ugene <- inferNetwork(Repressilator, mtry=3L)
test_that("cutoff too high", {
  expect_error(tuneThreshold(Repressilator, ugene, cutoffs = c(0.01, 0.99)))
})

# runs without errors
result <- tuneThreshold(Repressilator, ugene)
test_that("return values", {
  expect_equal(class(result), "ugene.analysis")
  expect_equal(length(result$stepErrors), 31)
  expect_equal(length(result$colErrors), 5)
  expect_equal(length(result$stepMasks), 31)
  expect_equal(length(result$colMasks), 5)
  expect_equal(dim(result$stepMasks[[1]])[1], 6)
  expect_equal(dim(result$stepMasks[[1]])[2], 6)
  expect_equal(dim(result$colMasks[[1]])[1], 6)
  expect_equal(dim(result$colMasks[[1]])[2], 6)
})

# [END]
