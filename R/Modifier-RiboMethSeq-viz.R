#' @include RNAmodR.RiboMethSeq.R
#' @include Modifier-RiboMethSeq-class.R
NULL

RNAMODR_RMS_PLOT_DATA <- c("ends",
                           "scoreA",
                           "scoreB",
                           "scoreRMS")
RNAMODR_RMS_PLOT_DATA_NAMES <- c(ends = "5'- & 3'-ends",
                                 scoreA = "Score A",
                                 scoreB = "Score B",
                                 scoreRMS = "Score RiboMethSeq")
RNAMODR_RMS_PLOT_DATA_COLOURS <- c(ends = "#FBB4AE",
                                   scoreA = "#B3CDE3",
                                   scoreB = "#CCEBC5",
                                   scoreRMS = "#DECBE4")
#' @rdname ModRiboMethSeq
#' 
#' @name visualizeData
#' @details 
#' \code{ModRiboMethSeq} specific arguments for \link{visualizeData}:
#' \itemize{
#' \item{\code{colour} - }{a named character vector of \code{length = 4} 
#' for the colours of the individual histograms. The names are expected to be 
#' \code{c("ends","scoreA","scoreB","scoreRMS")}}
#' }
NULL

#' @rdname ModRiboMethSeq
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModRiboMethSeq",
                        coord = "GRanges"),
  definition = function(x,
                        coord,
                        type = c("ends","scoreA","scoreB","scoreRMS"),
                        window.size = 15L,
                        ...) {
    if(missing(type)){
      type <- RNAMODR_RMS_PLOT_DATA
    }
    callNextMethod(x = x,
                   coord = coord,
                   type = type,
                   window.size = window.size,
                   ...)
  }
)
#' @rdname ModRiboMethSeq
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "ModRiboMethSeq"),
  definition = function(x,
                        name,
                        from,
                        to,
                        type = c("ends","scoreA","scoreB","scoreRMS"),
                        ...) {
    if(missing(type)){
      type <- RNAMODR_RMS_PLOT_DATA
    }
    callNextMethod(x = x,
                   name,
                   from,
                   to,
                   type = type,
                   ...)
  }
)

setMethod(
  f = ".dataTracks",
  signature = signature(x = "ModRiboMethSeq",
                        data = "GRanges",
                        seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x,
                        data,
                        seqdata,
                        sequence,
                        args) {
    requireNamespace("Gviz")
    n <- ncol(mcols(data))
    colour <- args[["colour"]]
    if(is.na(colour) || length(colour) != n){
      colour <- RNAMODR_RMS_PLOT_DATA_COLOURS
    }
    dts <- lapply(seq_len(n),
                  function(i){
                    column <- colnames(mcols(data)[i])
                    colour <- colour[column]
                    name <- RNAMODR_RMS_PLOT_DATA_NAMES[column]
                    dt <- DataTrack(data,
                              data = column,
                              name = name,
                              fill = colour,
                              type = "histogram")
                    if(column %in% c("scoreA","scoreRMS")){
                      displayPars(dt)$ylim = c(0,1)
                    }
                    displayPars(dt)$background.title <- "#FFFFFF"
                    displayPars(dt)$fontcolor.title <- "#000000"
                    displayPars(dt)$col.axis <- "#000000"
                    dt
                  })
    names(dts) <- colnames(mcols(data))
    dts
  }
)
