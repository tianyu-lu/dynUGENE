# test_that("Alpha inference not allowed", {
#   expect_error(inferSSNetwork(Repressilator, alpha = 'from.data'))
# })
#
# data <- grndata::syntren300.data
# ugene <- inferSSNetwork(data)
# test_that("Large example runs properly", {
#   expect_equal(dim(ugene$network)[1], 300)
#   expect_equal(dim(ugene$network)[2], 300)
#   expect_equal(sum(ugene$network > 1), 0)  # scores between 0 and 1
#   expect_equal(sum(ugene$network < 0), 0)
#   expect_equal(length(ugene$model), 300)
# })
#
# # [END]
