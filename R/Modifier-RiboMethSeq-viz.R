#' @include RNAmodR.RiboMethSeq.R
#' @include Modifier-RiboMethSeq-class.R
NULL

RNAMODR_RMS_PLOT_DATA <- c("ends",
                           "scoreA",
                           "scoreB",
                           "scoreRMS",
                           "scoreMean")
RNAMODR_RMS_PLOT_DATA_DEFAULT <- c("scoreRMS",
                                   "scoreMean")
RNAMODR_RMS_PLOT_DATA_NAMES <- c(ends = "5'- & 3'-ends",
                                 scoreA = "Score A",
                                 scoreB = "Score B",
                                 scoreRMS = "Score RiboMethSeq",
                                 scoreMean = "ScoreMean")
RNAMODR_RMS_PLOT_DATA_COLOURS <- c(ends = "#FBB4AE",
                                   scoreA = "#B3CDE3",
                                   scoreB = "#CCEBC5",
                                   scoreRMS = "#DECBE4",
                                   scoreMean = "#DECBE4")

.norm_viz_mod_rms_args <- function(input, type){
  if(!all(type %in% RNAMODR_RMS_PLOT_DATA)){
    stop("Type '",type,"' is not valid. Valid types are: '",
         paste0(RNAMODR_RMS_PLOT_DATA, collapse = "','"),"'.",
         call. = FALSE)
  }
  colour <- input[["colour"]]
  if(!is.null(input[["colour"]])){
    colour <- RNAmodR:::.norm_viz_colour(input[["colour"]], type)
  } else {
    colour <- RNAMODR_RMS_PLOT_DATA_COLOURS[type]
  }
  input <- list(type = type,
                colour = colour)
  input
}

#' @rdname ModRiboMethSeq-functions
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "ModRiboMethSeq"),
  definition = function(x, name, type, ...) {
    args <- .norm_viz_mod_rms_args(list(...), type)
    data <- RNAmodR:::.get_data_for_visualization(x, name)
    data <- unlist(data)
    dts <- lapply(
      args[["type"]],
      function(t){
        dt <- Gviz::DataTrack(range = data[,t],
                              groups = t,
                              name = RNAMODR_RMS_PLOT_DATA_NAMES[t],
                              col = args[["colour"]][t],
                              type = "histogram")
        if(t %in% c("scoreA","scoreRMS","scoreMean")){
          Gviz::displayPars(dt)$ylim = c(0,1)
        }
        Gviz::displayPars(dt)$background.title <- "#FFFFFF"
        Gviz::displayPars(dt)$fontcolor.title <- "#000000"
        Gviz::displayPars(dt)$col.axis <- "#000000"
        Gviz::displayPars(dt) <- args[names(args) != "type"]
        dt
      })
    names(dts) <- args[["type"]]
    dts
  }
)

#' @rdname ModRiboMethSeq-functions
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModRiboMethSeq",
                        coord = "GRanges"),
  definition = function(x, coord,
                        type = c("ends","scoreA","scoreB","scoreRMS","scoreMean"),
                        window.size = 15L,
                        ...) {
    if(missing(type)){
      type <- RNAMODR_RMS_PLOT_DATA_DEFAULT
    }
    type <- match.arg(type,c("scoreRMS","ends","scoreA","scoreB","scoreMean"),
                      several.ok = TRUE)
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
                        type = c("ends","scoreA","scoreB","scoreRMS","scoreMean"), ...) {
    if(missing(type)){
      type <- RNAMODR_RMS_PLOT_DATA_DEFAULT
    }
    type <- match.arg(type,c("scoreRMS","ends","scoreA","scoreB","scoreMean"),
                      several.ok = TRUE)
    callNextMethod(x = x, name = name, from = from, to = to, type = type, ...)
  }
)

#' @rdname ModRiboMethSeq-functions
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModSetRiboMethSeq",
                        coord = "GRanges"),
  definition = function(x, coord,
                        type = c("scoreRMS","ends","scoreA","scoreB","scoreMean"),
                        window.size = 15L,
                        ...) {
    if(missing(type)){
      type <- RNAMODR_RMS_PLOT_DATA_DEFAULT
    }
    type <- match.arg(type,c("scoreRMS","ends","scoreA","scoreB","scoreMean"),
                      several.ok = TRUE)
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
                        type = c("scoreRMS","ends","scoreA","scoreB","scoreMean"), ...) {
    if(missing(type)){
      type <- RNAMODR_RMS_PLOT_DATA_DEFAULT
    }
    type <- match.arg(type,c("scoreRMS","ends","scoreA","scoreB","scoreMean"),
                      several.ok = TRUE)
    callNextMethod(x = x, name = name, from = from, to = to, type = type, ...)
  }
)
