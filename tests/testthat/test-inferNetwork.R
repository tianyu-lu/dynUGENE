data <- read.csv("data/Repressilator.csv")
data$m1 <- NA

test_that("test errors", {
  expect_error( inferNetwork(data) )
})

data <- read.csv("data/Repressilator.csv")
ret <- inferNetwork(data)
test_that("return values", {
  expect_equal( class(ret), 'ugene' )
  expect_equal(sum(ret$network > 1), 0)  # scores between 0 and 1
  expect_equal(sum(ret$network < 0), 0)
  expect_equal(sum(ret$alpha > 1), 0)    # alpha (decay) can't be larger than 1
})
