#' @title RNAmodR.RiboMethSeq
#' 
#' @author Felix G M Ernst [aut], Denis L.J Lafontaine [fnd]
#' 
#' @description
#' `RNAmodR.RiboMethSeq` implements the detection of 2'-O methylations from 
#' RiboMethSeq data using the workflow and class the package `RNAmodR` provides.
#' 
#' @seealso Further details are described in the man pages of the 
#' \code{\link{Modifier}} object and the vignettes.
#'
#' @docType package
#' @name RNAmodR.RiboMethSeq
NULL

#' @import methods
#' @import RNAmodR
#' @import S4Vectors
#' @import BiocGenerics
#' @import IRanges
#' @import GenomicRanges
#' @import Gviz
NULL
requireNamespace("IRanges")
requireNamespace("RNAmodR")
