context("RiboMethSeq score A")
# .get_heatmap_data
test_that("RiboMethSeq score A:",{
  # score A
  # .calculate_ribometh_score_A
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_A(1,1,0,1,0)
  expect_equal(actual,0)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_A(100000,1,0,1,0)
  expect_equal(actual,0)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_A(1,10,0,10,0)
  expect_equal(actual,0.75)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_A(1,10,1,10,1)
  expect_equal(round(actual,7),0.7272727)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_A(1,100,0,100,0)
  expect_equal(round(actual,7),0.9705882)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_A(10,100,0,100,0)
  expect_equal(round(actual,7),0.8108108)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_A(10,100,10,100,10)
  expect_equal(round(actual,7),0.7920792)
})

context("RiboMethSeq score B")
# .get_heatmap_data
test_that("RiboMethSeq score B:",{
  # score B
  # .calculate_ribometh_score_B
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_B(1,1,0,1,0)
  expect_equal(actual,0)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_B(1,10,1,10,1)
  expect_equal(actual,4.5)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_B(1,100,1,100,1)
  expect_equal(actual,49.5)
  area <- c(100,50,200)
  weights <- c(1,1,1)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_B(1,area,weights,area,weights)
  expect_equal(round(actual,2),57.83)
  weights <- c(1,0.75,0.5)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_B(1,area,weights,area,weights)
  expect_equal(round(actual,2),52.28)
  weights <- c(1,0.5,0)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_B(1,area,weights,area,weights)
  expect_equal(round(actual,2),41.17)
})

context("RiboMethSeq score RMS")
# .get_heatmap_data
test_that("RiboMethSeq score RMS:",{
  # score RMS
  # .calculate_ribometh_score_meth
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_meth(1,1,0,1,0)
  expect_equal(actual,0)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_meth(1,10,1,10,1)
  expect_equal(actual,0.9)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_meth(1,100,1,100,1)
  expect_equal(actual,0.99)
  area <- c(100,50,200)
  weights <- c(1,1,1)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_meth(1,area,weights,area,weights)
  expect_equal(round(actual,7),0.9914286)
  weights <- c(1,0.75,0.5)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_meth(1,area,weights,area,weights)
  expect_equal(round(actual,7),0.9905263)
  weights <- c(1,0.5,0)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_meth(1,area,weights,area,weights)
  expect_equal(round(actual,7),0.988)
})


context("RiboMethSeq scoreMean")
# .get_heatmap_data
test_that("RiboMethSeq scoreMean:",{
  # score RMS
  # .calculate_ribometh_score_meth
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_mean(1,1)
  expect_equal(actual,0)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_mean(1,10)
  expect_equal(actual,0.9)
  actual <- RNAmodR.RiboMethSeq:::.calculate_ribometh_score_mean(1,100)
  expect_equal(actual,0.99)
})

# context("RiboMethSeq score MAX")
# # .get_heatmap_data
# test_that("RiboMethSeq score MAX:",{
#             # score MAX
#             # .calculate_ribometh_score_max
#           })