#' @include RiboMethSeq.R
NULL

#' @name ModRiboMethSeq
#' @aliases RiboMethSeq ModRiboMethSeq
#' 
#' @title ModRiboMethSeq
#' @description 
#' title
#' 

NULL

#' @rdname ModRiboMethSeq
#' @export
setClass("ModRiboMethSeq",
         contains = c("Modifier"),
         prototype = list(mod = "I",
                          dataType = "PileupSequenceData"))


setMethod(
  f = "initialize", 
  signature = signature(.Object = "ModRiboMethSeq"),
  definition = function(.Object,
                        bamfiles,
                        fasta,
                        gff) {
    .Object <- callNextMethod(.Object,
                              bamfiles,
                              fasta,
                              gff)
    return(.Object)
  }
)

# constructors -----------------------------------------------------------------


setGeneric( 
  name = "ModRiboMethSeq",
  def = function(x,
                 ...) standardGeneric("ModRiboMethSeq")
)
# Create Modifier class from file character, fasta and gff file
#' @rdname ModRiboMethSeq
#' @export
setMethod("ModRiboMethSeq",
          signature = c(x = "character"),
          function(x,
                   fasta,
                   gff,
                   modifications = NULL,
                   ...){
            x
          }
)

# Create Modifier class from bamfiles, fasta and gff file
#' @rdname ModRiboMethSeq
#' @export
setMethod("ModRiboMethSeq",
          signature = c(x = "BamFileList"),
          function(x,
                   fasta,
                   gff,
                   modifications = NULL,
                   ...){
            x
          }
)

# Create Modifier class from existing SequenceData
#' @rdname ModRiboMethSeq
#' @export
setMethod("ModRiboMethSeq",
          signature = c(x = "SequenceData"),
          function(x,
                   modifications = NULL,
                   ...){
            x
          }
)

# functions --------------------------------------------------------------------

.get_rms_scores <- function(data){
  
}

.find_rms <- function(mod,
                          data,
                          args){
}

#' @rdname ModRiboMethSeq
#' @export
setMethod("modify",
          signature = c(x = "ModRiboMethSeq"),
          function(x,
                   ...){
            x
          }
)
