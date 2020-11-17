rates <- estimateDecayRates(Repressilator)
test_that("decay rates between 0 and 1", {
  expect_equal(sum(rates > 1), 0)
  expect_equal(sum(rates < 0), 0)
})

# other tests are integrated in test-inferNetwork.R and test-inferSSNetwork.R
