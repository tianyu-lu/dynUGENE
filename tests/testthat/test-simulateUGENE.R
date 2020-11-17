ugene <- inferNetwork(Repressilator, mtry=3L)
x0 <- Repressilator[1, 2:7]
x0.flat <- unlist(x0)

test_that("input correctness", {
  expect_error(simulateUGENE(Repressilator, x0))
  expect_error(simulateUGENE(ugene, x0.flat))
  expect_error(simulateUGENE(ugene, x0, mask = c(0.1, 0.2)))
})

trajectory <- simulateUGENE(ugene, x0)
test_that("output correctness", {
  expect_equal(class(trajectory), "simulation")
  expect_equal(length(trajectory$t), 1001)
  expect_equal(dim(trajectory$x)[1], 1001)
  expect_equal(dim(trajectory$x)[2], 6)
})

test_that("integration with plotTrajectory", {
  expect_silent(plotTrajectory(trajectory, c("p3", "p2", "p1")))
})
