#' @title Imputing missing values
#'
#' @description
#' Imputes missing values or values of low quality.
#'
#' @param ms An S4-object of class \link{BCmeasurement}
#' @param method Imputation method. See details.
#' @param max.missing Portion of values per metabolite that are in maximum
#' allowed to be missing in order to be included in imputation.
#' @param keep.status Status of values that should be kept as learning data (i.e., not
#' subject to imputation).
#'
#' @returns An object of class \link{BCmeasurement} with imputed values. If a
#' value was imputed or not is stored in the `object@status` slot with the tag
#' "imputed".
#'
#' @details
#' Possible imputation methods are:
#' * <b>LOD</b> - Values with the status "< LOD" are imputed with the respective LOD
#' of the focal metabolite and Biocrates plate.
#' * `RF` - Nonparametric Missing Value Imputation using Random Forest as
#' implemented in \link[missForest]{missForest}. Note that RF imputation in
#' non-deterministic.
#'
#' @export
imputeBC <- function(ms, method, max.missing = 0.2,
                     keep.status = c("Valid","< LLOQ","> ULOQ", "Imputed")) {

  # Construct sub-table, where all NAs will be imputed
  dat <- ms@conc
  dat[which(!(ms@status %in% keep.status), arr.ind = TRUE)] <- NA_real_
  metind_keep <- which(apply(dat,2,function(x) sum(is.na(x))/length(x)) <= max.missing)
  dat <- dat[, metind_keep]
  ind_imp <- which(is.na(dat), arr.ind = TRUE)

  # Random Forest Imputation
  if(method == "RF") {
    dat <- imputeBC_RF(dat)
  }

  # update status
  ind_imp2 <- ind_imp; ind_imp2[,2] <- metind_keep[ind_imp2[,2]]
  ms@status[ind_imp2] <- "Imputed"

  # Update concentration table
  ms@conc[ind_imp2] <- dat[ind_imp]

  return(ms)
}


#' @import missForest
imputeBC_RF <- function(dat) {
  dati <- missForest(dat, maxiter = 100)
  return(dati$ximp)
}

