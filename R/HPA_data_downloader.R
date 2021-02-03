#' @title HPA_data_downloader
#'
#' @description Used to download required data for HPAStainR
#'
#' @param tissue_type A character string that determines which HPA data you
#'  want to download from the website. Has to be both' (default), 'normal', 
#'  or 'cancer'.
#' @param save_file A boolean determining if you want the HPA data downloaded
#'  permanently or temporarily. Default is TRUE, meaning the file will be saved
#'   in the given 'save_location', default being the current working directory.
#' @param save_location A character string indicating where you want the files
#'  to be saved if you are saving them. If the file(s) already exists in that
#'  location, those will be loaded instead of redownloading the files.
#' @param version_date_normal a character string indicating what date of the 
#' normal tissue file you want returned if you have saved one. If you say "last"
#' it will return the most recently downloaded file, otherwise give a date in
#' YYYY-MM-DD format.
#' @param version_date_cancer a character string indicating what date of the 
#' cancer tissue file you want returned if you have saved one. If you say "last"
#' it will return the most recently downloaded file, otherwise give a date in
#' YYYY-MM-DD format.
#' 
#' @return  List of dataframes or dataframe depending on tissue_type arguement.
#'  If tissue_type == 'both' it will be a list of dataframes.
#'
#' @examples
#' HPA_data <- HPA_data_downloader(tissue_type = 'both', save_file = FALSE)
#' ## Access normal data
#' HPA_data$hpa_dat
#' ## Access cancer data
#' HPA_data$cancer_dat
#'
#'
#' ## Download only the normal tissue data
#' HPA_normal_data <- HPA_data_downloader('normal', save_file = FALSE)
#' @importFrom utils download.file read.table unzip
#' @importFrom data.table fread
#' @export
HPA_data_downloader <- function(tissue_type = c("both", "normal", "cancer"), 
    save_file = TRUE, save_location = "", version_date_normal = "last", version_date_cancer = "last") {
    
    ## A suggest function from Bioconductor allowing for better resource
    ## querying
    getURL <- function(URL, FUN, ..., N.TRIES = 1L) {
        N.TRIES <- as.integer(N.TRIES)
        stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))
        
        while (N.TRIES > 0L) {
            result <- tryCatch(FUN(URL, ...), error = identity)
            if (!inherits(result, "error")) 
                break
            N.TRIES <- N.TRIES - 1L
        }
        
        if (N.TRIES == 0L) {
            stop("'getURL()' failed:", "\n  URL: ", URL, "\n  error: ", 
                conditionMessage(result))
        }
        
        result
    }
    
    tissue_type <- tissue_type[1]

    ## normal URL
    norm_url <- "https://www.proteinatlas.org/download/normal_tissue.tsv.zip"
    ## path URL
    path_url <- "https://www.proteinatlas.org/download/pathology.tsv.zip"
    
        
    ## First section runs if files are already downloaded to save time
    
    ## Get date 
    cur_date <- Sys.Date()
    
    
    ## Make if else here that searches other files in the area and picks the one
    ## with the most recent date.
    dir_files <- ifelse(save_location == "", ".", save_location)
    
    if (version_date_normal == "last") {
       
        ## Find most recent normal file
        norm_files <- list.files(dir_files)[grepl("normal_tissue",
                                                  list.files(dir_files))]
        version_date_normal <- max(str_sub(norm_files, 15, 24))
    }
    if (version_date_cancer == "last") {
        
        ## Find most recent cancer file
        canc_files <- list.files(dir_files)[grepl("pathology",
                                                  list.files(dir_files))]
        version_date_cancer <- max(str_sub(norm_files, 11, 20))
    }
    
    ## Saving strings for normal tissue and cancer to keep consistency
    ## throughout the code
    saved_normal_dat <- paste0(save_location,
                         "normal_tissue_", version_date_normal,".tsv.zip")
    saved_cancer_dat <- paste0(save_location, 
                               "pathology_", version_date_cancer,".tsv.zip")
    
    if ((tissue_type == "both" | tissue_type == "normal") &
        file.exists(saved_normal_dat)) {
        ## Normal tissue
        hpa_dat <- data.table::fread(
            unzip(saved_normal_dat), 
            header = TRUE, sep = "\t", stringsAsFactors = FALSE,
            data.table = FALSE, check.names = TRUE)
        
        if (tissue_type != "both") {
            return(hpa_dat)
        }
    }
    
    if ((tissue_type == "both" | tissue_type == "cancer") &
        file.exists(saved_cancer_dat)) {
        ## Cancer tissue
        
        cancer_dat <- data.table::fread(
            unzip(saved_cancer_dat), 
            header = TRUE, sep = "\t", stringsAsFactors = FALSE,
            data.table = FALSE, check.names = TRUE)
        
        if (tissue_type != "both") {
            return(cancer_dat)
        }
    }
    if (tissue_type == "both" & file.exists(
        saved_cancer_dat) & 
        file.exists(saved_normal_dat)) {
        all_dat <- list()
        all_dat$hpa_dat <- hpa_dat
        all_dat$cancer_dat <- cancer_dat
        return(all_dat)
    }
    
    
    
    ## Begins section where files are downloaded 
    
    if (save_file == FALSE) {
        if (tissue_type == "both" | tissue_type == "normal") {
            ## Normal tissue
            ## 2/3/2021 Switched to fread require tmp_dir and unzip code
            ## as fread can't read files made using unz()
            temp <- tempfile()
            tmp_dir = tempdir()
            getURL(URL = norm_url, 
                FUN = download.file, destfile = temp)
            unzip(temp,  "normal_tissue.tsv", exdir = tmp_dir)
            hpa_dat <- data.table::fread(file.path(tmp_dir, "normal_tissue.tsv"),
                                  header = TRUE,
                                  sep = "\t", stringsAsFactors = FALSE,
                                  data.table = FALSE, check.names = TRUE)
            unlink(temp)
            if (tissue_type != "both") {
                return(hpa_dat)
            }
        }
        
        if (tissue_type == "both" | tissue_type == "cancer") {
            ## Cancer tissue
            temp <- tempfile()
            tmp_dir = tempdir()
            getURL(URL = path_url, 
                FUN = download.file, destfile = temp)
            unzip(temp,  "pathology.tsv", exdir = tmp_dir)
            cancer_dat <- data.table::fread(file.path(tmp_dir, "pathology.tsv"),
                                            header = TRUE, 
                sep = "\t", stringsAsFactors = FALSE,
                data.table = FALSE, check.names = TRUE)
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
            
            getURL(URL = norm_url, 
                FUN = download.file, destfile = paste0(save_location,
                                                       "normal_tissue.tsv.zip"))
            hpa_dat <- data.table::fread(unzip(paste0(save_location,
                                               "normal_tissue.tsv.zip")), 
                header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                data.table = FALSE, check.names = TRUE)
            
            if (tissue_type != "both") {
                return(hpa_dat)
            }
        }
        
        if (tissue_type == "both" | tissue_type == "cancer") {
            ## Cancer tissue
            
            getURL(URL = path_url, 
                FUN = download.file,
                destfile = paste0(save_location, "pathology.tsv.zip"))
            cancer_dat <- data.table::fread(unzip(paste0(save_location,
                                                  "pathology.tsv.zip")), 
                header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                data.table = FALSE, check.names = TRUE)
            
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
