#' @include RNAmodR.RiboMethSeq.R
NULL

#' @name ModRiboMethSeq
#' @aliases RiboMethSeq ModRiboMethSeq ModSetRiboMethSeq
#' 
#' @author Felix G.M. Ernst
#' 
#' @title ModRiboMethSeq
#' @description 
#' Among the various post-transcriptional RNA modifications, 2'-O methylations
#' are quite common in rRNA and tRNA. They confere resistance to alkaline 
#' degradation by preventing a nucleophilic attack on the 3'-phosphate 
#' especially in flexible RNA, which is fascilitated by high pH conditions.
#' This property can be queried using a method called RiboMethSeq (Birkedahl et 
#' al. 2015, Marchand et al. 2017) for which RNA is treated in alkaline 
#' conditions and RNA fragments are used to prepare a sequencing library.
#' 
#' At position containing a 2'-O methylations, read ends are less frequent, 
#' which is used to detect and score the2'-O methylations.
#' 
#' The \code{ModRiboMethSeq} class uses the the 
#' \code{\link[RNAmodR:ProtectedEndSequenceData]{ProtectedEndSequenceData}}
#' class to store and aggregate data along the transcripts. The calculated 
#' scores follow the nomenclature of Birkedahl et al. (2015) with the names
#' \code{scoreRMS} (default), \code{scoreA} and \code{scoreB}.
#' 
#' The score MAX as described by Marchand et al. (2017) are not implemented, 
#' yet, since an unambigeous description is not available from the literature.
#' 
#' @param x the input which can be of the different types depending on whether
#' a \code{ModRiboMethSeq} or a \code{ModSetRiboMethSeq} object is to be 
#' constructed. For more information have a look at the documentation of
#' the \code{\link[RNAmodR:Modifier]{Modifier}} and 
#' \code{\link[RNAmodR:ModifierSet]{ModifierSet}} classes.
#' @param annotation annotation data, which must match the information contained
#' in the BAM files. This is parameter is only required if \code{x} if not a 
#' \code{Modifier} object.
#' @param sequences sequences matching the target sequences the reads were 
#' mapped onto. This must match the information contained in the BAM files. This
#' is parameter is only required if \code{x} if not a \code{Modifier} object.
#' @param seqinfo An optional \code{\link[GenomeInfoDb:Seqinfo-class]{Seqinfo}} 
#' argument or character vector, which can be coerced to one, to subset the 
#' sequences to be analyzed on a per chromosome basis.
#' @param ... Optional arguments overwriting default values, which are
#' \itemize{
#' \item{weights:} {The weights used for calculating the scores B and RMS 
#' (default: \code{weights = c(0.9,1,0,1,0.9)})}
#' \item{flankingRegion:} {The size of the flanking region used for calculation 
#' of score A as an integer value (default: \code{flankingRegion = 6L})}
#' \item{minSignal:} {The minimal singal at the position as integer value 
#' (default: \code{minSignal = 10L}). If the reaction is very specific a lower
#' value and even 0L may need to be used}
#' \item{minScoreA:} {minimum for score A to identify 2'-O methylated positions 
#' de novo (default: \code{minScoreA = 0.6})}
#' \item{minScoreB:} {minimum for score B to identify 2'-O methylated positions 
#' de novo (default: \code{minScoreB = 3.0})}
#' \item{minScoreRMS:} {minimum for score RMS to identify 2'-O methylated 
#' positions de novo (default: \code{minScoreRMS = 0.75})}
#' \item{scoreOperator:} {how the minimal score should be used as logical 
#' operator. "&" requires all minimal values to be exceeded, whereas "|" detects
#' positions, if at least one minimal values is exceeded (default: 
#' \code{scoreOperator = "|"})}
#' \item{maxLength:} {The default read length. Reads with this length or longer
#' are discarded, since they represent non-fragemented reads. This might need to
#' be adjusted for individual samples dending on the experimental conditions.
#' This is argument is passed on to 
#' \code{\link[RNAmodR:ProtectedEndSequenceData]{ProtectedEndSequenceData}}
#' (default: \code{maxLength = 50L})}
#' \item{other arguments}{which are passed on to 
#' \code{\link[RNAmodR:ProtectedEndSequenceData]{ProtectedEndSequenceData}}}
#' }
#' To disable minimal values for modification calling, set them to \code{0}.
#' It is not advised to set them all to \code{0}.
#' 
#' @references 
#' - Birkedal U, Christensen-Dalsgaard M, Krogh N, Sabarinathan R, Gorodkin J, 
#' Nielsen H (2015): "Profiling of ribose methylations in RNA by high-throughput
#' sequencing." Angewandte Chemie (International ed. in English) 54 (2), 
#' P. 451–455. DOI: 
#' \href{https://doi.org/10.1002/anie.201408362}{10.1002/anie.201408362}.
#' 
#' - Marchand V, Ayadi L, El Hajj A, Blanloeil-Oillo F, Helm M, Motorin Y 
#' (2017): "High-Throughput Mapping of 2'-O-Me Residues in RNA Using 
#' Next-Generation Sequencing (Illumina RiboMethSeq Protocol)." Methods in 
#' molecular biology (Clifton, N.J.) 1562, P. 171–187. DOI: 
#' \href{https://doi.org/10.1007/978-1-4939-6807-7_12}{10.1007/978-1-4939-6807-7_12}.
NULL

#' @rdname ModRiboMethSeq
#' @export
setClass("ModRiboMethSeq",
         contains = c("Modifier"),
         prototype = list(mod = c("Am","Cm","Gm","Um"),
                          score = "scoreRMS",
                          dataType = "ProtectedEndSequenceData"))

# settings ---------------------------------------------------------------------

#' @name ModRiboMethSeq-functions
#' @aliases aggregate modify settings visualizeData visualizeDataByCoord
#' 
#' @title Functions for ModRiboMethSeq
#' 
#' @description
#' All of the functions of \code{\link[RNAmodR:Modifier]{Modifier}} and the
#' \code{\link[RNAmodR:Modifier]{ModifierSet}} classes are inherited by the 
#' \code{ModRiboMethSeq} and \code{ModSetRiboMethSeq} classes.
#' 
#' @param x a \code{\link[RNAmodR:Modifier]{Modifier}} or a 
#' \code{\link[RNAmodR:ModifierSet]{ModifierSet}} object. For more details see 
#' also the man pages for the functions mentioned below.
#' @param value See \code{\link[RNAmodR:Modifier]{settings}}
#' @param force See \code{\link[RNAmodR:aggregate]{aggregate}}
#' @param coord,name,from,to,type,window.size,... See 
#' \code{\link[RNAmodR:visualizeData]{visualizeData}}
#' 
#' @details 
#' \code{ModRiboMethSeq} specific arguments for \link{visualizeData}:
#' \itemize{
#' \item{\code{colour} - }{a named character vector of \code{length = 4} 
#' for the colours of the individual histograms. The names are expected to be 
#' \code{c("ends","scoreA","scoreB","scoreRMS")}}
#' }
#' 
#' @importMethodsFrom RNAmodR modify aggregate settings visualizeData 
#' visualizeDataByCoord
NULL

.valid_rms_weights <- function(weights){
  if(!is.numeric(weights) | !is.atomic(weights)){
    return(FALSE)
  }
  if((length(weights) %% 2) != 1){
    return(FALSE)
  }
  if(length(weights) < 3L){
    return(FALSE)
  }
  if(weights[ceiling(length(weights)/2)] != 0){
    return(FALSE)
  }
  TRUE
}

.norm_rms_args <- function(input){
  maxLength <- 50L # for all scores
  minSignal <- 10L # for all scores
  flankingRegion <- 6L # for score A
  minScoreA <- 0.6 # for score A
  minScoreB <- 4.0 # for score B
  minScoreRMS <- 0.75 # for score C/RMS
  weights <- c(0.9,1,0,1,0.9) # for score B/C/RMS
  scoreOperator <- "&"
  if(!is.null(input[["weights"]])){
    weights <- input[["weights"]]
    if(!.valid_rms_weights(weights)){
      stop("'weights' must be a numeric vector of uneven length. The middle ",
           "position will be used as the current position.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["minSignal"]])){
    minSignal <- input[["minSignal"]]
    if(!is.integer(minSignal) | minSignal < 1L){
      stop("'minSignal' must be integer with a value higher than 0.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["minScoreA"]])){
    minScoreA <- input[["minScoreA"]]
    if(!is.numeric(minScoreA) | minScoreA < 0 | minScoreA > 1){
      stop("'minScoreA' must be numeric with a value between 0 and 1.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["flankingRegion"]])){
    flankingRegion <- input[["flankingRegion"]]
    if(!is.integer(flankingRegion) | flankingRegion < 1L){
      stop("'flankingRegion' must be integer with a value higher than 0L.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["minScoreB"]])){
    minScoreB <- input[["minScoreB"]]
    if(!is.numeric(minScoreB) | minScoreB < 0){
      stop("'minScoreB' must be numeric and greater then 0.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["minScoreRMS"]])){
    minScoreRMS <- input[["minScoreRMS"]]
    if(!is.numeric(minScoreRMS) | minScoreRMS < 0 | minScoreRMS > 1){
      stop("'minScoreRMS' must be numeric with a value between 0 and 1.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["scoreOperator"]])){
    scoreOperator <- input[["scoreOperator"]]
    if(!(scoreOperator %in% c("|","&"))){
      stop("'scoreOperator' must be either '|' or '&'.",
           call. = FALSE)
    }
  }
  args <- RNAmodR:::.norm_args(input)
  args <- c(args,
            list(maxLength = maxLength,
                 weights = weights,
                 minSignal = minSignal,
                 flankingRegion = flankingRegion,
                 minScoreA = minScoreA,
                 minScoreB = minScoreB,
                 minScoreRMS = minScoreRMS,
                 scoreOperator = scoreOperator))
  args
}

#' @rdname ModRiboMethSeq-functions
#' @export
setReplaceMethod(f = "settings", 
                 signature = signature(x = "ModRiboMethSeq"),
                 definition = function(x, value){
                   x <- callNextMethod()
                   value <- .norm_rms_args(value)
                   x@arguments[names(value)] <- unname(value)
                   x
                 })


# constructor ------------------------------------------------------------------

#' @rdname ModRiboMethSeq
#' @export
ModRiboMethSeq <- function(x, annotation = NA, sequences = NA, seqinfo = NA,
                           ...){
  RNAmodR:::Modifier("ModRiboMethSeq", x = x, annotation = annotation,
                     sequences = sequences, seqinfo = seqinfo, ...)
}


# functions --------------------------------------------------------------------

# calculates score A according to Birkedal et al. 2015
# it is simplified to use the mean/sd for the neighboring positions
# and not the left and right mean/sd seperatly
.calculate_ribometh_score_A_c <- function(n,
                                          meanL,
                                          sdL,
                                          meanR,
                                          sdR){
  dividend <- (2 * n  + 1)
  divisor <- (0.5 * abs(meanL - sdL)) + n + (0.5 * abs(meanR - sdR)) + 1
  ans <- 1 - (dividend / divisor)
  ans <- vapply(ans,max,numeric(1),0)
  ans[is.na(ans)] <- 0
  return(ans)
}
.calculate_ribometh_score_A <- compiler::cmpfun(.calculate_ribometh_score_A_c)

# calculates score B according to Birkedal et al. 2015
# it is simplified to use the weighted neighboring positions
# and not the left and right area seperatly
.calculate_ribometh_score_B_c <- function(n,
                                          areaL,
                                          weightsL,
                                          areaR,
                                          weightsR){
  waL <- sum(weightsL * areaL) / sum(weightsL)
  waR <- sum(weightsR * areaR) / sum(weightsR)
  dividend <- abs(n - 0.5 * ( waL + waR ) )
  divisor <- (n + 1)
  ans <- dividend / divisor
  ans[is.na(ans)] <- 0
  return(ans)
}
.calculate_ribometh_score_B <- compiler::cmpfun(.calculate_ribometh_score_B_c)

# calculates score C according to Birkedal et al. 2015
# it is simplified to use the weighted neighboring positions
# and not the left and right area seperatly
.calculate_ribometh_score_meth_c <- function(n,
                                             areaL,
                                             weightsL,
                                             areaR,
                                             weightsR){
  waL <- sum(weightsL * areaL) / sum(weightsL)
  waR <- sum(weightsR * areaR) / sum(weightsR)
  dividend <- n
  divisor <- 0.5 * ( waL + waR )
  ans <- 1 - (dividend / divisor)
  ans <- vapply(ans,max,numeric(1),0)
  ans[is.na(ans)] <- 0
  return(ans)
}
.calculate_ribometh_score_meth <- compiler::cmpfun(.calculate_ribometh_score_meth_c)

# calculates score MAX according to Marchand et al. 2016
# not used since no clear description available in the literature
# .calculate_ribometh_score_max_c <- function(n,
#                                              area,
#                                              weights){
#   browser()
# }
# .calculate_ribometh_score_max <- compiler::cmpfun(.calculate_ribometh_score_max_c)


# RiboMeth scores --------------------------------------------------------------

# calculates score A according to Birkedal et al. 2015
.get_score_A <- function(data,
                         countsL,
                         countsR){
  n <- seq_along(data)
  countsMeanL <- IRanges::NumericList(
    lapply(n,
           function(j){
             unlist(lapply(countsL[[j]],
                           mean,
                           na.rm = TRUE))
           }))
  countsSdL <- IRanges::NumericList(
    lapply(n,
           function(j){
             unlist(lapply(countsL[[j]],
                           BiocGenerics::sd,
                           na.rm = TRUE))
           }))
  countsSdL[is.na(countsSdL)] <- 0
  countsMeanR <- IRanges::NumericList(
    lapply(n,
           function(j){
             unlist(lapply(countsR[[j]],
                           mean,
                           na.rm = TRUE))
           }))
  countsSdR <- IRanges::NumericList(
    lapply(n,
           function(j){
             unlist(lapply(countsR[[j]],
                           BiocGenerics::sd,
                           na.rm = TRUE))
           }))
  countsSdR[is.na(countsSdR)] <- 0
  # calc score per replicate
  scoreA <- IRanges::NumericList(mapply(FUN = .calculate_ribometh_score_A,
                                        data,
                                        countsMeanL,
                                        countsSdL,
                                        countsMeanR,
                                        countsSdR))
  return(scoreA)
}

# calculates score B according to Birkedal et al. 2015
.get_score_B <- function(data,
                         countsL,
                         weightsL,
                         countsR,
                         weightsR){
  scoreB <- IRanges::NumericList(mapply(FUN = .calculate_ribometh_score_B,
                                        data,
                                        countsL,
                                        weightsL,
                                        countsR,
                                        weightsR))
  return(scoreB)
}

# calculates score C according to Birkedal et al. 2015
.get_score_meth <- function(data,
                            countsL,
                            weightsL,
                            countsR,
                            weightsR){
  scoreMeth <- IRanges::NumericList(mapply(FUN = .calculate_ribometh_score_meth,
                                           data,
                                           countsL,
                                           weightsL,
                                           countsR,
                                           weightsR))
  return(scoreMeth)
}

# calculates score C according to Marchand et al. 
# .get_score_max <- function(data,
#                            countsL,
#                            weightsL,
#                            countsR,
#                            weightsR){
#   browser()
#   scoreMAX <- IRanges::NumericList(mapply(FUN = .calculate_ribometh_score_max,
#                                  data,
#                                  countsL,
#                                  weightsL,
#                                  countsR,
#                                  weightsR))
#   return(scoreMAX)
# }

.norm_rms_weights <- function(x){
  weights <- settings(x,"weights")
  if(is.null(names(weights))){
    names(weights) <- seq_along(weights) - ((length(weights) / 2) + 0.5)
  }
  names <- as.integer(names(weights))
  if(any(is.na(names))){
    stop("If 'weights' is named, all names must be coercible to integer ",
         "values.",
         call. = FALSE)
  }
  if(length(which(names == 0L)) != 1){
    stop("If 'weights' is named, exactly one name must be '0'",
         call. = FALSE)
  }
  weights
}

.aggregate_rms <- function(x){
  message("Aggregating data and calculating scores...")
  # parameter data
  weights <- .norm_rms_weights(x)
  weightPositions <- as.integer(names(weights))
  # ToOo check for continuity
  if(length(weights) != length(weightPositions)){
    stop("Something went wrong.")
  }
  # get the means. the sds arecurrently disregarded for this analysis
  mod <- aggregate(seqData(x), condition = "Treated")
  means <- IRanges::IntegerList(mod@unlistData[,which(grepl("mean",
                                                            colnames(mod@unlistData)))])
  means@partitioning <- mod@partitioning
  # set up variables
  n <- length(mod)
  nV <- seq_len(n)
  lengths <- S4Vectors::lengths(mod)
  pos <- lapply(lengths,seq_len)
  flankingRegion <- settings(x,"flankingRegion")
  positionsR <- seq_len(flankingRegion)
  positionsL <- rev(positionsR) * -1
  weightPositionsL <- weightPositions[weightPositions < 0]
  weightPositionsR <- weightPositions[weightPositions > 0]
  weightPositionsC <- which(weightPositions == 0)
  # subset to neightbouring positions based on the size of flankingRegions
  neighborCountsLFR <- 
    lapply(nV,
           function(j){
             IRanges::IntegerList(lapply(pos[[j]],
                                         function(k){
                                           f <- positionsL + k
                                           ans <- means[[j]][f[f > 0]]
                                           ans <- ans[!is.na(ans)]
                                           ans
                                         }))
           })
  neighborCountsRFR <- 
    lapply(nV,
           function(j){
             IRanges::IntegerList(lapply(pos[[j]],
                                         function(k){
                                           f <- positionsR + k
                                           ans <- means[[j]][f[f > 0]]
                                           ans <- ans[!is.na(ans)]
                                           ans
                                         }))
           })
  # subset to neightbouring positions based on the size of the weights
  neighborCountsL <- 
    lapply(nV,
           function(j){
             IRanges::IntegerList(lapply(pos[[j]],
                                         function(k){
                                           f <- weightPositionsL + k
                                           ans <- means[[j]][f[f > 0]]
                                           ans <- ans[!is.na(ans)]
                                           ans
                                         }))
           })
  neighborCountsR <- 
    lapply(nV,
           function(j){
             IRanges::IntegerList(lapply(pos[[j]],
                                         function(k){
                                           f <- weightPositionsR + k
                                           ans <- means[[j]][f[f > 0]]
                                           ans <- ans[!is.na(ans)]
                                           ans
                                         }))
           })
  # create list of weights vector alongside the neighbor counts
  weightsListL <- 
    lapply(nV,
           function(j){
             IRanges::NumericList(mapply(
               function(k,l){
                 f <- weightPositionsL + k
                 ans <- weights[which(f > 0 & f <= l)]
                 ans
               },
               pos[[j]],
               MoreArgs = list(l = lengths[j])))
           })
  weightsListR <- 
    lapply(nV,
           function(j){
             IRanges::NumericList(mapply(
               function(k,l,of){
                 f <- weightPositionsR + k
                 f <- which(f > 0 & f <= l) + of
                 ans <- weights[f]
                 ans
               },
               pos[[j]],
               MoreArgs = list(l = lengths[j],
                               of = weightPositionsC)))
           })
  # calculate the actual scores
  scoreA <- .get_score_A(means,
                         neighborCountsLFR,
                         neighborCountsRFR)
  scoreB <- .get_score_B(means,
                         neighborCountsL,
                         weightsListL,
                         neighborCountsR,
                         weightsListR)
  scoreRMS <- .get_score_meth(means,
                              neighborCountsL,
                              weightsListL,
                              neighborCountsR,
                              weightsListR)
  # scoreMAX <- .get_score_max()
  ans <- S4Vectors::DataFrame(ends = unlist(means),
                              scoreA = unlist(scoreA),
                              scoreB = unlist(scoreB),
                              scoreRMS = unlist(scoreRMS),
                              row.names = NULL)
  ans <- IRanges::SplitDataFrameList(ans)
  ans@partitioning <- mod@partitioning
  ans
}

#' @rdname ModRiboMethSeq-functions
#' @export
setMethod(
  f = "aggregate", 
  signature = signature(x = "ModRiboMethSeq"),
  definition = 
    function(x, force = FALSE){
      if(missing(force)){
        force <- FALSE
      }
      if(!hasAggregateData(x) || force){
        x@aggregate <- .aggregate_rms(x)
        x@aggregateValidForCurrentArguments <- TRUE
      }
      x
    }
)

.get_rms_scores <- function(data){
  list(score = data$scoreRMS,
       scoreA = data$scoreA,
       scoreB = data$scoreB)
}

.find_rms <- function(x){
  message("Searching for 2'-O methylations...")
  #
  letters <- IRanges::CharacterList(strsplit(as.character(sequences(x)),""))
  grl <- ranges(x)
  # get the aggregate data
  mod <- aggregateData(x)
  # setup args
  minSignal <- settings(x,"minSignal")
  minScoreA <- settings(x,"minScoreA")
  minScoreB <- settings(x,"minScoreB")
  minScoreRMS <- settings(x,"minScoreRMS")
  scoreOperator <- settings(x,"scoreOperator")
  # find modifications
  modifications <- mapply(
    function(m,l,r){
      rownames(m) <- seq_len(BiocGenerics::width(r))
      m <- m[!is.na(m$scoreA) &
               !is.na(m$scoreB) &
               !is.na(m$scoreRMS),]
      if(nrow(m) == 0L) return(NULL)
      m <- m[m$ends >= minSignal,]
      if(nrow(m) == 0L) return(NULL)
      m <- m[mapply(Reduce,
                    rep(scoreOperator,nrow(m)),
                    m$scoreA >= minScoreA,
                    m$scoreB >= minScoreB,
                    m$scoreRMS >= minScoreRMS),]
      if(nrow(m) == 0L) return(NULL)
      ans <- RNAmodR::.constructModRanges(
        r,
        m,
        modType = "Am",
        scoreFun = RNAmodR.RiboMethSeq:::.get_rms_scores,
        source = "RNAmodR.RiboMethSeq",
        type = "RNAMOD")
      ans$mod <- paste0(l[BiocGenerics::start(ans)],"m")
      ans
    },
    mod,
    letters,
    grl)
  modifications <- GenomicRanges::GRangesList(
    modifications[!vapply(modifications,
                          is.null,
                          logical(1))])
  unname(unlist(modifications))
}

#' @rdname ModRiboMethSeq-functions
#' @export
setMethod("modify",
          signature = c(x = "ModRiboMethSeq"),
          function(x, force){
            # get the aggregate data
            x <- aggregate(x, force)
            x@modifications <- .find_rms(x)
            x@modificationsValidForCurrentArguments <- TRUE
            message("done.")
            x
          }
)


# ModSetRiboMethSeq ------------------------------------------------------------

#' @rdname ModRiboMethSeq
#' @export
setClass("ModSetRiboMethSeq",
         contains = "ModifierSet",
         prototype = list(elementType = "ModRiboMethSeq"))

#' @rdname ModRiboMethSeq
#' @export
ModSetRiboMethSeq <- function(x, annotation = NA, sequences = NA, seqinfo = NA){
  RNAmodR::ModifierSet("ModRiboMethSeq", x, annotation = annotation,
                       sequences = sequences, seqinfo = seqinfo)
}
