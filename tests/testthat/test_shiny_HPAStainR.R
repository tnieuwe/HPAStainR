context("Shiny App Test")
library(HPAStainR)
module <- shiny_HPAStainR(hpa_dat = HPA_data_downloader(save_file = FALSE)[[1]],
                          cancer_dat = HPA_data_downloader(save_file = FALSE)[[2]])

test_that("Is this a shiny object?", {
  expect_is(module, "shiny.appobj")
})