#' Structure of the S4 class "BCmeasurement"
#'
#' @aliases BCmeasurement
#'
#' @exportClass BCmeasurement
#'
#' @slot mets data.table with column information (i.e., metabolite info)
#' @slot samples data.table with row information (i.e., sample info)
#' @slot conc Numeric matrix of measured concentrations. Rows are samples,
#' columns are metabolites.
#' @slot indicators Numeric matrix with concentration-derived metabolism
#' indicator metrics. (Not yet implemented)
#' @slot status Character matrix of each metabolites measurement quality.

setClass("BCmeasurement",

         slots = c(
           mets = "data.table",
           samples = "data.table",
           conc = "matrix",
           indicators = "matrix", # TODO
           status = "matrix"
         ))

setMethod("initialize", "BCmeasurement",
          function(.Object,
                   mets,
                   samples,
                   conc,
                   indicators,
                   status,
                   ...) {
            .Object <- callNextMethod(.Object, ...)

            .Object@mets <- mets
            .Object@samples <- samples
            .Object@conc <- conc
            .Object@indicators <- indicators
            .Object@status <- status

            return(.Object)
          }
)
