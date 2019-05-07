context("RiboMethSeq")
# .get_heatmap_data
test_that("RiboMethSeq:",{
  # weight validation
  expect_false(RNAmodR.RiboMethSeq:::.valid_rms_weights(c(1,1,1)))
  expect_true(RNAmodR.RiboMethSeq:::.valid_rms_weights(c(1,0,1)))
  expect_false(RNAmodR.RiboMethSeq:::.valid_rms_weights(c(1,1,0,1)))
  expect_false(RNAmodR.RiboMethSeq:::.valid_rms_weights(c(1,0)))
  expect_false(RNAmodR.RiboMethSeq:::.valid_rms_weights(c(1)))
  expect_true(RNAmodR.RiboMethSeq:::.valid_rms_weights(c(1,1,0,1,1)))
  # argument normalization
  input <- list()
  actual <- RNAmodR.RiboMethSeq:::.norm_rms_args(input)
  expect_type(actual,"list")
  expect_named(actual,c("minCoverage",
                        "minReplicate",
                        "find.mod",
                        "maxLength",
                        "weights",
                        "minSignal",
                        "flankingRegion",
                        "minScoreA",
                        "minScoreB",
                        "minScoreRMS",
                        "minScoreMean",
                        "flankingRegionMean",
                        "scoreOperator"))
  expect_error(RNAmodR.RiboMethSeq:::.norm_rms_args(list(weights = c(1,1,0,1))))
  expect_error(RNAmodR.RiboMethSeq:::.norm_rms_args(list(minSignal = 10)))
  expect_equal(RNAmodR.RiboMethSeq:::.norm_rms_args(list(minSignal = 10L)),
               actual)
  expect_error(RNAmodR.RiboMethSeq:::.norm_rms_args(list(minScoreA = 2L)))
  expect_equal(RNAmodR.RiboMethSeq:::.norm_rms_args(list(minScoreA = 0.6)),
               actual)
  expect_error(RNAmodR.RiboMethSeq:::.norm_rms_args(list(flankingRegion = 2.5)))
  expect_equal(RNAmodR.RiboMethSeq:::.norm_rms_args(list(flankingRegion = 6L)),
               actual)
  expect_error(RNAmodR.RiboMethSeq:::.norm_rms_args(list(minScoreB = "a")))
  expect_equal(RNAmodR.RiboMethSeq:::.norm_rms_args(list(minScoreB = 3.6)),
               actual)
  expect_error(RNAmodR.RiboMethSeq:::.norm_rms_args(list(minScoreRMS = 2L)))
  expect_equal(RNAmodR.RiboMethSeq:::.norm_rms_args(list(minScoreRMS = 0.75)),
               actual)
  expect_error(RNAmodR.RiboMethSeq:::.norm_rms_args(list(scoreOperator = 2L)))
  expect_equal(RNAmodR.RiboMethSeq:::.norm_rms_args(list(scoreOperator = "&")),
               actual)
  
  data <- data.frame(ends = "1", scoreRMS = "a", scoreA = "b", scoreB = "c",
                     scoreMean = "d",
                     stringsAsFactors = FALSE)
  actual <- RNAmodR.RiboMethSeq:::.get_rms_scores(data)
  expect_type(actual,"list")
  expect_named(actual,c("score","scoreA","scoreB","scoreMean"))
})
