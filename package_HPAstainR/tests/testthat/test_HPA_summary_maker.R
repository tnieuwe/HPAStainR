context("HPA Summary maker")
library(HPAStainR)
test_that("HPA summary make makes a tibble",{
          expect_is(HPA_summary_maker(
            HPA_data_downloader(save_file = FALSE)[[1]]),
            "tbl_df")
  })
