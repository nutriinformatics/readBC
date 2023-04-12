#' @title Reads an excel file with a Biocrates metabolomics data set
#'
#' @description
#' A function that transforms a excel table with Biocrates results into an
#' S4 R object that is easier to handle for data processing, QC, imputation,
#' statistics and data merging.
#'
#' @param file Path to excel file.
#' @param sheet Name of data sheet with main table.
#' @param beautify logical. Should specific metabolite names and classes be
#' renamed? See details. (Default: TRUE)
#'
#' @details
#' If TRUE, the parameter 'beautify' renames the following entities:
#' * "Amine Oxides" -> "Amine oxides"
#' * "Aminoacids" -> "Amino acids"
#' * "Aminoacids Related" -> "Amino acid-related"
#' * "Bile Acids" -> "Bile acids"
#' * "Biogenic Amines" -> "Biogenic amines"
#' * "Carboxylic Acids" -> "Carboxylic acids"
#' * "Cholesterol Esters" -> "Cholesterol esters"
#' * "Fatty Acids" -> "Fatty acids"
#' * "Indoles Derivatives" -> "Indole derivatives"
#' * "Nucleobases Related" -> "Nucleobase-related"
#' * "Vitamins & Cofactors" -> "Vitamins and cofactors"
#'
#' @returns An object of class \link{BCmeasurement}.
#'
#' @examples
#' metdat <- readBC(system.file("extdata", "DZHK_Serum_Feb2021_readBC.xlsx",
#'                              package="readBC"),
#'                  sheet = "Data Export")
#'
#' @import tidyxl
#' @import data.table
#'
#' @export
readBC <- function(file, sheet = "Data Export", beautify = TRUE) {

  # read all cell info from excel file
  sheet_name <- sheet
  cells <- xlsx_cells(file)
  cells <- data.table(cells)
  cells <- cells[sheet == sheet_name]

  #---------#
  # margins #
  #---------#

  # concentrations
  conc_s <- c(cells[character == "TMAO",row] + 1,
              cells[character == "Class",col] + 1)
  tmp <- cells[character == "Class",row]
  conc_e <- c(cells[col == 10, max(row)],
              cells[row == tmp & character != "Metabolism Indicators",max(col)])
  rm(tmp)

  # col info (metabolites)
  mets_s <- cells[character == "Class",c(row,col)]
  mets_e <- c(conc_s[1]-1,
              conc_e[2])

  # row info (samples)
  samples_s <- c(conc_s[1],
                 cells[row == conc_s[1] - 1, min(col)])
  samples_e <- c(conc_e[1],
                 conc_s[2]-1)

  #--------------#
  # extract data #
  #--------------#

  # col info (metabolites)
  mets <- data.table(Class = cells[row == mets_s[1] & col > mets_s[2] & col <= mets_e[2], character])
  for(i in (mets_s[1]+1):mets_e[1]) {
    cname <- cells[row == i & col == mets_s[2], character]
    if(i == mets_e[1])
      cname <- "Metabolite"

    if(grepl("^LOD|^ULOQ|^LLOQ", cname)) {
      mets[[cname]] <- rep(NA_real_, nrow(mets))
      tmp_ind <- cells[row == i & col > mets_s[2] & col <= mets_e[2], col] - mets_s[2]
      mets[[cname]][tmp_ind] <- cells[row == i & col > mets_s[2] & col <= mets_e[2], numeric]
      rm(tmp_ind)
    } else {
      mets[[cname]] <- rep(NA_character_, nrow(mets))
      tmp_ind <- cells[row == i & col > mets_s[2] & col <= mets_e[2], col] - mets_s[2]
      mets[[cname]][tmp_ind] <- cells[row == i & col > mets_s[2] & col <= mets_e[2], character]
    }
  }

  # row info (samples)
  samples <- data.table(tmp = samples_s[1]:samples_e[1])
  for(j in samples_s[2]:samples_e[2]) {
    cname <- cells[row == samples_s[1]-1 & col == j, character]
    samples[[cname]] <- ifelse(cells[row >= samples_s[1] & row <= samples_e[1] & col == j, is.na(numeric)],
                               cells[row >= samples_s[1] & row <= samples_e[1] & col == j, character],
                               cells[row >= samples_s[1] & row <= samples_e[1] & col == j, numeric])
  }
  samples[, tmp := NULL]

  # measured concentrations
  conc <- cells[row >= conc_s[1] & col >= conc_s[2] & row <= conc_e[1] & col <= conc_e[2], numeric]
  conc <- matrix(conc, byrow = TRUE,
                 nrow = conc_e[1] - conc_s[1] + 1,
                 ncol = conc_e[2] - conc_s[2] + 1)

  # measurement status
  status <- cells[row >= conc_s[1] & col >= conc_s[2] & row <= conc_e[1] & col <= conc_e[2], character]
  status <- matrix(status, byrow = TRUE,
                   nrow = conc_e[1] - conc_s[1] + 1,
                   ncol = conc_e[2] - conc_s[2] + 1)

  OPs <- colnames(mets); OPs <- OPs[grep("^OP ",OPs)]
  plate_codes <- unique(samples$`Plate Bar Code`)

  for(j in 1:ncol(conc)){
    cat("\r",j,"/", ncol(conc))

    vPlates <- t(mets[j, ..OPs]); vPlates <- rownames(vPlates)[!is.na(vPlates[,1])]
    vPlates <- gsub("^OP ","", vPlates)

    for(k in plate_codes) {
      plateIDs  <- unlist(strsplit(k, " \\| "))
      plateIDs <- intersect(plateIDs, vPlates)
      plateRegEx <- paste(c(plateIDs, gsub("-","/",plateIDs)), collapse = "|")
      ind_spls <- which(samples$`Plate Bar Code` == k)

      relMetInfo <- colnames(mets)[grepl(plateRegEx, colnames(mets))]

      # LLOQ
      tmp <- relMetInfo[grep("^LLOQ",relMetInfo)]
      LLOQ <- as.numeric(mets[j, ..tmp]); LLOQ <- LLOQ[!is.na(LLOQ)]

      # ULOQ
      tmp <- relMetInfo[grep("^ULOQ",relMetInfo)]
      ULOQ <- as.numeric(mets[j, ..tmp]); ULOQ <- ULOQ[!is.na(ULOQ)]
      if(ULOQ == 0)
        ULOQ <- Inf

      status[ind_spls,j] <- ifelse(is.na(conc[ind_spls,j]),
                                   status[ind_spls,j],
                                   ifelse(conc[ind_spls,j] < LLOQ,
                                          "< LLOQ",
                                          ifelse(conc[ind_spls,j] > ULOQ,
                                                 "> ULOQ",
                                                 ifelse(!is.na(status[ind_spls,j]),
                                                        status[ind_spls,j],
                                                        "Valid"))))

    }
  }

  #----------#
  # Beautify #
  #----------#
  if(beautify) {
    mets[Class == "Amine Oxides", Class := "Amine oxides"]
    mets[Class == "Aminoacids", Class := "Amino acids"]
    mets[Class == "Aminoacids Related", Class := "Amino acid-related"]
    mets[Class == "Bile Acids", Class := "Bile acids"]
    mets[Class == "Biogenic Amines", Class := "Biogenic amines"]
    mets[Class == "Carboxylic Acids", Class := "Carboxylic acids"]
    mets[Class == "Cholesterol Esters", Class := "Cholesterol esters"]
    mets[Class == "Fatty Acids", Class := "Fatty acids"]
    mets[Class == "Indoles Derivatives", Class := "Indole derivatives"]
    mets[Class == "Nucleobases Related", Class := "Nucleobase-related"]
    mets[Class == "Vitamins & Cofactors", Class := "Vitamins and cofactors"]
  }

  #------------#
  # Finalizing #
  #------------#
  res <- new("BCmeasurement",
             mets, samples, conc, matrix(NA_real_, nrow = nrow(conc), ncol = 0),
             status)

}
