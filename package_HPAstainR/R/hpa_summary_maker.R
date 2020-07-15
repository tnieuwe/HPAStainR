#' @title hpa_summary_maker
#'
#' @description Used to generate a summary file used in the second tab of the Shiny app version of HPAStainR
#'
#' @param tissue_type A character string that determines which HPA data you want to download from the website. Has to be
#' "both" (default), "normal", or "cancer".
#' @param save_file A boolean determining if you want the HPA data downloaded permanently or temporarily. Default is TRUE,
#' meaning the file will be saved in the given "save_location", default being the current working directory.
#' @param save_location A character string indicating where you want the files to be saved if you are saving them.
#' If the file(s) already exists in that location, those will be loaded instead of re-downloading the files.
#'
#' @return  A dataframe summarizing the amount of proteins tested to detected, used for the shiny app.
#'
#' @examples
#' ## Generate the summarized HPA file
#' hpa_summary <- hpa_summary_maker(hpa_dat)
#' @export
#'
#'


hpa_summary_maker <- function(hpa_dat) {
  hpa_dat_new <- hpa_dat %>%
    mutate(
      Tissue = gsub("[[:digit:]]+", "", Tissue),
      tissue_cell = paste0(toupper(str_trim(Tissue)), " - ", Cell.type)
    ) %>%
    distinct() %>%
    group_by(tissue_cell) %>%
    summarise(proteins = n(), detected = sum(!(Level %in% "Not detected"))) %>%
    mutate(det_o_test = detected / proteins)
}
