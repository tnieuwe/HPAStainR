#' @title HPA_summary_maker
#'
#' @description Used to generate a summary file used in the second tab of the Shiny app version of HPAStainR
#'
#' @param hpa_dat The dataframe of normal tissue data downloaded by HPA_data_downloader()
#'
#' @return  A dataframe summarizing the amount of proteins tested to detected, used for the shiny app.
#'
#' @examples
#' ## Load in data from downloader
#' HPA_data <- HPA_data_downloader(save_file = FALSE)
#' ## Generate the summarized HPA file
#' hpa_summary <- HPA_summary_maker(HPA_data$hpa_dat)
#' @importFrom dplyr distinct n
#' @export
#'
#'


HPA_summary_maker <- function(hpa_dat) {
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
