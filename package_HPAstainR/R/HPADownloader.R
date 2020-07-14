#' @title HPA_data_downloader
#'
#' @description Used to download required data for HPAStainR
#'
#' @param tissue_type A character string that determines which HPA data you want to download from the website. Has to be
#' "both" (default), "normal", or "cancer".
#' @param save_file A boolean determining if you want the HPA data downloaded permanently or temporarily. Default is TRUE,
#' meaning the file will be saved in the given "save_location", default being the current working directory.
#' @param save_location A character string indicating where you want the files to be saved if you are saving them.
#' If the file(s) already exists in that location, those will be loaded instead of redownloading the files.
#'
#' @return  List of dataframes or dataframe depending on tissue_type arguement. If tissue_type == "both" it will be a list of dataframes.
#'
#' @examples
#' HPA_data <- HPA_data_downloader()
#' ## Access normal data
#' HPA_data$hpa_dat
#' ## Access cancer data
#' HPA_data$cancer_dat
#'
#'
#' ## Download only the normal tissue data
#' HPA_normal_data <- HPA_data_downloader("normal")
#' @export
HPA_data_downloader <- function(tissue_type = c("both", "normal", "cancer"),
                                save_file = TRUE,
                                save_location = "") {
  tissue_type <- tissue_type[1]

  ## First section runs if files are already downloaed to save time
  if ((tissue_type == "both" | tissue_type == "normal") & file.exists(paste0(save_location, "normal_tissue.tsv.zip"))) {
    ## Normal tissue
    hpa_dat <- read.table(unzip(paste0(save_location, "normal_tissue.tsv.zip")),
      header = TRUE, sep = "\t", stringsAsFactors = FALSE
    )

    if (tissue_type != "both") {
      return(hpa_dat)
    }
  }

  if ((tissue_type == "both" | tissue_type == "cancer") & file.exists(paste0(save_location, "pathology.tsv.zip"))) {
    ## Cancer tissue

    cancer_dat <- read.table(unzip(paste0(save_location, "pathology.tsv.zip")),
      header = TRUE, sep = "\t", stringsAsFactors = FALSE
    )

    if (tissue_type != "both") {
      return(cancer_dat)
    }
  }
  if (tissue_type == "both" & file.exists(paste0(save_location, "pathology.tsv.zip")) &
    file.exists(paste0(save_location, "normal_tissue.tsv.zip"))) {
    all_dat <- list()
    all_dat$hpa_dat <- hpa_dat
    all_dat$cancer_dat <- cancer_dat
    return(all_dat)
  }

  if (save_file == F) {
    if (tissue_type == "both" | tissue_type == "normal") {
      ## Normal tissue
      temp <- tempfile()
      download.file("https://www.proteinatlas.org/download/normal_tissue.tsv.zip", temp)
      hpa_dat <- read.table(unz(temp, "normal_tissue.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      unlink(temp)
      if (tissue_type != "both") {
        return(hpa_dat)
      }
    }

    if (tissue_type == "both" | tissue_type == "cancer") {
      ## Cancer tissue
      temp <- tempfile()
      download.file("https://www.proteinatlas.org/download/pathology.tsv.zip", temp)
      cancer_dat <- read.table(unz(temp, "pathology.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      unlink(temp)
      if (tissue_type != "both") {
        return(cancer_dat)
      }
    }
    if (tissue_type == "both") {
      all_dat <- list()
      all_dat$hpa_dat <- hpa_dat
      all_dat$cancer_dat <- cancer_dat
      return(all_dat)
    }
  }
  if (save_file == TRUE) {
    if (tissue_type == "both" | tissue_type == "normal") {
      ## Normal tissue

      download.file("https://www.proteinatlas.org/download/normal_tissue.tsv.zip",
        destfile = paste0(save_location, "normal_tissue.tsv.zip")
      )
      hpa_dat <- read.table(unzip(paste0(save_location, "normal_tissue.tsv.zip")),
        header = TRUE, sep = "\t", stringsAsFactors = FALSE
      )

      if (tissue_type != "both") {
        return(hpa_dat)
      }
    }

    if (tissue_type == "both" | tissue_type == "cancer") {
      ## Cancer tissue

      download.file("https://www.proteinatlas.org/download/pathology.tsv.zip",
        destfile = paste0(save_location, "pathology.tsv.zip")
      )
      cancer_dat <- read.table(unzip(paste0(save_location, "pathology.tsv.zip")),
        header = TRUE, sep = "\t", stringsAsFactors = FALSE
      )

      if (tissue_type != "both") {
        return(cancer_dat)
      }
    }
    if (tissue_type == "both") {
      all_dat <- list()
      all_dat$hpa_dat <- hpa_dat
      all_dat$cancer_dat <- cancer_dat
      return(all_dat)
    }
  }
}
