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

#' @rdname ModRiboMethSeq-functions
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModRiboMethSeq",
                        coord = "GRanges"),
  definition = function(x, coord,
                        type = c("ends","scoreA","scoreB","scoreRMS"),
                        window.size = 15L,
                        ...) {
    if(missing(type)){
      type <- RNAMODR_RMS_PLOT_DATA
    }
    callNextMethod(x = x, coord = coord, type = type, window.size = window.size,
                   ...)
  }
)
#' @rdname ModRiboMethSeq-functions
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "ModRiboMethSeq"),
  definition = function(x, name, from, to,
                        type = c("ends","scoreA","scoreB","scoreRMS"), ...) {
    if(missing(type)){
      type <- RNAMODR_RMS_PLOT_DATA
    }
    callNextMethod(x = x, name, from, to, type = type, ...)
  }
)

setMethod(
  f = ".dataTracks",
  signature = signature(x = "ModRiboMethSeq",
                        data = "GRanges",
                        seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x, data, seqdata, sequence, args) {
    n <- ncol(S4Vectors::mcols(data))
    colour <- args[["colour"]]
    if(is.na(colour) || length(colour) != n){
      colour <- RNAMODR_RMS_PLOT_DATA_COLOURS[seq.int(1,n)]
    }
    if(is.null(names(colour))){
      names(colour) <- colnames(S4Vectors::mcols(data))
    }
    dts <- lapply(seq_len(n),
                  function(i){
                    column <- colnames(S4Vectors::mcols(data)[i])
                    colour <- colour[column]
                    name <- RNAMODR_RMS_PLOT_DATA_NAMES[column]
                    dt <- Gviz::DataTrack(data,
                                          data = column,
                                          name = name,
                                          fill = colour,
                                          type = "histogram")
                    if(column %in% c("scoreA","scoreRMS")){
                      Gviz::displayPars(dt)$ylim = c(0,1)
                    } else if(!is.null(args[["ylim"]])){
                      Gviz::displayPars(dt)$ylim <- args[["ylim"]]
                    }
                    Gviz::displayPars(dt)$background.title <- "#FFFFFF"
                    Gviz::displayPars(dt)$fontcolor.title <- "#000000"
                    Gviz::displayPars(dt)$col.axis <- "#000000"
                    if(!is.null(args[["data.track.pars"]])){
                      Gviz::displayPars(dt) <- args[["data.track.pars"]]
                    }
                    dt
                  })
    names(dts) <- colnames(S4Vectors::mcols(data))
    dts
  }
)

#' @rdname ModRiboMethSeq-functions
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModSetRiboMethSeq",
                        coord = "GRanges"),
  definition = function(x, coord,
                        type = c("scoreRMS","ends","scoreA","scoreB"),
                        window.size = 15L,
                        ...) {
    if(missing(type)){
      type <- "scoreRMS"
    }
    type <- match.arg(type,c("scoreRMS","ends","scoreA","scoreB"))
    callNextMethod(x = x, coord = coord, type = type, window.size = window.size,
                   ...)
  }
)
#' @rdname ModRiboMethSeq-functions
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "ModSetRiboMethSeq"),
  definition = function(x, name, from, to,
                        type = c("scoreRMS","ends","scoreA","scoreB"), ...) {
    if(missing(type)){
      type <- "scoreRMS"
    }
    type <- match.arg(type,c("scoreRMS","ends","scoreA","scoreB"))
    callNextMethod(x = x, name, from, to, type = type, ...)
  }
)
