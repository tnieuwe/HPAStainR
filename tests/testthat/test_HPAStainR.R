context("HPAStainR")
library(HPAStainR)
test_that("Does HPAStainR run?",{
  expect_is(HPAStainR(
    gene_list = c("PRSS1", "CELA3A", "PNLIP", "PRL"),
    hpa_dat = HPA_data_downloader(save_file = FALSE)[[1]],
    cancer_dat = HPA_data_downloader(save_file = FALSE)[[2]]),
    "tbl_df")
      })
test_that("Columns of the top row are not NAs", {
  expect_equal(sum(is.na((HPAStainR(
    gene_list = c("PRSS1", "CELA3A", "PNLIP", "PRL"),
    hpa_dat = HPA_data_downloader(save_file = FALSE)[[1]],
    cancer_dat = HPA_data_downloader(save_file = FALSE)[[2]]))[1,])),
    0)
})

