#' @title RNAmodR.RiboMethSeq
#' 
#' @author Felix G M Ernst [aut], Denis L J Lafontaine [fnd]
#' 
#' @description
#' `RNAmodR.RiboMethSeq` implements the detection of 2'-O methylations from 
#' RiboMethSeq data using the workflow and class the package `RNAmodR` provides.
#' 
#' @seealso Further details are described in the man pages of the 
#' \code{\link[RNAmodR:Modifier-class]{Modifier}} object and the vignettes.
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
requireNamespace("RNAmodR")

#' @name RNAmodR.RiboMethSeq-datasets
#' @title Example data in the RNAmodR.RiboMethSeq package
#' @description 
#' This contains an example ModifierSet object of type ModSetRiboMethSeq
#' @docType data
#' @usage msrms
#' @format a \code{ModSetRiboMethSeq} instance
#' @keywords datasets
"msrms"