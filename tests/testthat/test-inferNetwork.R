data <- Repressilator
data$m1 <- NA

test_that("does not accept data with NA", {
  expect_error( inferNetwork(data) )
})

ret <- inferNetwork(Repressilator)
test_that("check return values", {
  expect_equal( class(ret), 'ugene' )
  expect_equal(dim(ret$network)[1], 6) # 6x6 matrix (6 components in data)
  expect_equal(dim(ret$network)[2], 6)
  expect_equal(sum(ret$network > 1), 0)  # scores between 0 and 1
  expect_equal(sum(ret$network < 0), 0)
  expect_equal(sum(ret$alpha > 1), 0)    # alpha (decay) can't be larger than 1
  expect_equal(length(ret$model), 6) # six components in data = six random forests
  expect_equal(class(ret$model[[1]]), "randomForest")
})

# mask 3:
#   V1 V2 V3 V4 V5 V6
# 1 NA NA NA NA NA NA
# 2 NA NA  1  1 NA NA
# 3 NA NA  1 NA NA NA
# 4 NA NA NA NA  1  1
# 5 NA NA NA NA  1 NA
# 6  1  1 NA NA NA NA

mask3 <- as.matrix(read.csv("mask3.csv", row.names = 1))
mask3.invalid <- mask3
mask3.invalid[2,3] <- 0.89
ret <- inferNetwork(Repressilator, mask=mask3)
x0 <- Repressilator[1, 2:7]
test_that("masking works", {
  # using an input size of 6 should give an error when training data has < 6 nodes
  expect_error(predict(ret$model[[1]], x0))
  # masks can only contain 1 or NA
  expect_error(inferNetwork(Repressilator, mask=mask3.invalid))
})

test_that("multiple experiments error", {
  # Repressilator dataset does not have multiple experiments
  expect_error(inferNetwork(Repressilator, multiple.exp = TRUE))
})

# [END]
