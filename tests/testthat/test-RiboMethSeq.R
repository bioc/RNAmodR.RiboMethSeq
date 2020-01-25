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
                        "minSignal",
                        "flankingRegion",
                        "minScoreA",
                        "minScoreB",
                        "minScoreRMS",
                        "minScoreMean",
                        "flankingRegionMean",
                        "weights",
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
  # settings
  data("msrms", package = "RNAmodR.RiboMethSeq")
  expect_type(settings(msrms),"list")
  expect_type(settings(msrms[[1]]),"list")
  expect_error(settings(msrms[[1]]) <- 1,
               "'value' has to be a named.")
  expect_error(settings(msrms[[1]]) <- c(1),
               "'value' has to be a named.")
  # settings(msrms[[1]]) <- c(minCoverage = 11L)
  # expect_equal(settings(msrms[[1]])$minCoverage, 11L)
  # # aggregate Data
  expect_s4_class(sequenceData(msrms[[1]]), "ProtectedEndSequenceData")
  mod <- aggregate(sequenceData(msrms[[1]]), condition = "Treated")
  expect_named(mod)
  expect_s4_class(mod,"CompressedSplitDFrameList")
  expect_equal(colnames(mod[[1L]]),c("means.treated","sds.treated"))
  expect_equal(mod[[1L]][1,1], 14152)
  #
  actual <- RNAmodR.RiboMethSeq:::.aggregate_rms(msrms[[1]]) 
  expect_named(actual)
  expect_s4_class(actual,"CompressedSplitDFrameList")
  expect_equivalent(actual[[1]][1,3], 0)
  expect_equivalent(actual[[1]][2,3], 0.75416916223702946)
  expect_equal(colnames(actual[[1]]),c("ends","scoreA","scoreB","scoreRMS",
                                       "scoreMean"))
  # findMod
  actual <- findMod(msrms[[1]])
  expect_s4_class(actual,"GRanges")
  expect_length(actual,0)
  expect_equal(dim(mcols(actual)),c(0,0))
  actual <- findMod(msrms[[2]])
  expect_s4_class(actual,"GRanges")
  expect_length(actual,1)
  expect_equal(colnames(mcols(actual)),c("mod","source","type","score",
                                         "scoreA","scoreB","scoreMean",
                                         "Parent"))
  expect_equal(unique(mcols(actual)$mod),c("Um"))
  # show
  # expect_output(expect_warning(show(msrms),
  # "Settings were changed after data aggregation or"))
  # expect_warning(modifications(msrms),
  # "Settings were changed after data aggregation or")
  expect_true(all(modifications(msrms)[[1]] == actual))
  # plotting
  expect_error(RNAmodR.RiboMethSeq:::.norm_viz_mod_rms_args(list(), "1"),
               "Type '1' is not valid")
  actual <- RNAmodR.RiboMethSeq:::.norm_viz_mod_rms_args(list(), "scoreA")
  expect_type(actual,"list")
  expect_named(actual,c("type","colour"))
  expect_equal(actual[["type"]],"scoreA")
  actual <- RNAmodR.RiboMethSeq:::.norm_viz_mod_rms_args(list(), "scoreB")
  expect_equal(actual[["type"]],"scoreB")
  actual <- RNAmodR.RiboMethSeq:::.norm_viz_mod_rms_args(list(), "ends")
  expect_equal(actual[["type"]],"ends")
  # getDataTrack
  expect_error(getDataTrack(msrms[[1]]),
               'argument "type" is missing, with no default')
  expect_error(getDataTrack(msrms[[1]],type="score"),
               "Type 'score' is not valid")
  expect_error(getDataTrack(msrms[[1]],type="scoreB"),
               'argument "name" is missing')
  actual <- getDataTrack(msrms[[1]],name="1",type="scoreB")
  expect_type(actual,"list")
  expect_named(actual,"scoreB")
  expect_s4_class(actual[[1]],"DataTrack")
  # plotData
  actual <- plotData(msrms[[1]],name="1",type="scoreB")
  expect_type(actual,"list")
  expect_named(actual,c("Score B","RNASequenceTrack","titles"))
  expect_s4_class(actual[[1]],"DataTrack")
  expect_s4_class(actual[[2]],"RNASequenceTrack")
  expect_s4_class(actual[[3]],"ImageMap")
  actual <- plotData(msrms,name="1",type="scoreB")
  expect_type(actual,"list")
  expect_s4_class(actual[[1]],"DataTrack")
  expect_s4_class(actual[[2]],"DataTrack")
  expect_s4_class(actual[[3]],"RNASequenceTrack")
  expect_s4_class(actual[[4]],"ImageMap")
})
