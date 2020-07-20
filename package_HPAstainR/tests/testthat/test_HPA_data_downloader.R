context("HPA Downloader")
library(HPAStainR)

test_that("Files can be downloaded from the website", {
  expect_is(HPA_data_downloader(save_file = FALSE), "list")
  expect_is(HPA_data_downloader(save_file = FALSE)[[1]], "data.frame")
  expect_is(HPA_data_downloader(save_file = FALSE)[[2]], "data.frame")
})

test_that("File columns haven't changed", {
  expect_equal(colnames(HPA_data_downloader(save_file = FALSE)[[1]]),
               c("Gene",
                 "Gene.name",
                 "Tissue",
                 "Cell.type",
                 "Level",
                 "Reliability")
               )
  
  expect_equal(colnames(HPA_data_downloader(save_file = FALSE)[[2]]),
               c("Gene",
                 "Gene.name",
                 "Cancer",
                 "High",
                 "Medium",
                 "Low",
                 "Not.detected",
                 "prognostic...favourable",
                 "unprognostic...favourable",
                 "prognostic...unfavourable", 
                 "unprognostic...unfavourable")
  )
})
