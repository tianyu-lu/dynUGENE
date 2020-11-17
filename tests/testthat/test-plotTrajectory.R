x0 <- Repressilator[1, 2:7]
ugene <- inferNetwork(Repressilator)
trajectory <- simulate(ugene, x0)

test_that("input correctness", {
  expect_error(plotTrajectory(trajectory, c("I don't exist", "me neither")))
  expect_error(plotTrajectory(trajectory, c("m1","m1","m1","m1","m1","m1","m1")))
  expect_error(plotTrajectory(trajectory, c("t")))
})
