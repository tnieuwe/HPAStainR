#' @title shiny_HPAStainR
#'
#' @description Runs HPAStainR as a the shiny app found at
#'  https://32tim32.shinyapps.io/HPAStainR/
#'
#' @param hpa_dat A required dataframe that has the normal tissue dataframe
#'  (see HPA_data_downloader).
#' 
#' @param cancer_dat A required dataframe that has the cancer tissue dataframe
#'  (see HPA_data_downloader).
#'
#' @param cell_type_data An optional dataframe that comes out of the
#'  hpa_summary_maker function, only needed if you want the second tab of
#'  HPAStainR, which shows the ratio of tested proteins to stained proteins,
#'  to be functional.
#'
#'
#' @return  A locally ran shiny app
#'
#' @examples
#' ## Load in data from downloader
#' HPA_data <- HPA_data_downloader(save_file = FALSE)
#' ## Generate the summarized HPA file
#' hpa_summary <- HPA_summary_maker(HPA_data$hpa_dat)
#' ## Run with summary, commented out so example doesn't run indefinitely
#' ## shiny_HPAStainR(HPA_data$hpa_dat, HPA_data$cancer_dat, hpa_summary)
#' @importFrom utils write.csv
#' @export


shiny_HPAStainR <- function(hpa_dat, cancer_dat, cell_type_data = NULL) {
  shinyApp(
    ui = fluidPage(
      titlePanel("HPAStain.R"),
      # Option selection
      sidebarLayout(
        # List of genes
        sidebarPanel(
          paste("Provide a list of proteins (genes), of standard human",
          "nomenclature, and HPAStain.R will output the cell type they are",
          "most associated with based on patterns of staining in the Human",
          "Protein Atlas. Check the second tab for information on how often",
          "tissues are tested for proteins"),
          textAreaInput("gene_list",
                        paste("List of Proteins (comma separated,",
                        "space separated, or line separated):"),
                        "PRSS1,PNLIP,CELA3A,PRL"),
          submitButton(text = "Run HPAStain.R"),
          downloadButton("downloadData",
                         "Download Output as CSV"),
          selectInput("tissue_level",
                      "Report tissue source of cell type?",
                      c("Yes", "No")),
          selectInput("cancer_analysis",
                      "Normal/Cancer Tissue", c("normal", "both", "cancer")),
          selectInput("percent_or_count",
                      "Report counts, percents or both for expression levels?",
                      c("percent", "count", "both")),
          ## This scales the genes fo the percentage so you don't get 100
          ## enrichment scores for tissues where 1 protein was tested
          checkboxInput("scale_abundance",
                        "Scale results for proteins that have available data",
                        TRUE),
          selectInput("stringency",
                      "Confidence level of HPA data",
                      c("Normal", "High", "Low"),
                      selected = "Normal"),
          checkboxInput("tested_protein_column",
                        "Have a column of tested proteins",
                        TRUE),
          checkboxInput("stain_gene_results",
                        "Have a column of detected proteins",
                        TRUE),
          checkboxInput("adjusted_pvals",
                        "Include adjusted P-values",
                        TRUE),
          #
          numericInput("round_to",
                       "Round values to:",
                       value = 2,
                       min = 2,
                       max = 5),
          checkboxInput("csv_names",
                        "CSV friendly column name format",
                        FALSE),
          checkboxInput("drop_na_rows_in",
                        "Remove rows with no staining data",
                        FALSE),

          strong(textOutput("gene_list")),
          p(
            paste(
            "    HPAStain.R is an R based tool used to query the Human Protein",
            "Atlas for staining data."),
            paste(
            "The purpose of this tool is to test if a list of proteins is",
            "associated with a certain cell type in a tissue."),
            paste(
            "E.g. you have a list of protein coding genes from a differential",
            "expression single cell analysis and want to see if these proteins",
            "are associated with a known cell type."),
            paste(
            "Instead of querying HPA multiple times you can load your list in",
            "HPAStain.R which will return a ranked table of the cell types",
            "with the most protein staining."),
          
            ),
          br(),
          p(
            paste(
            "    HPAStain.R is limited to the data which is available on the",
            "HPA website, as a result, not all tissues or cell types are",
            "stained or characterized for your proteins of interest are in the",
            "database."),
            paste(
            "When data is lacking on a gene it will not be included in the",
            "column `tested proteins`, and when tissue data is missing it is",
            "simply left blank."
            )
          ),
          br(),
          "The data used is from the Human Protein Atlas.",
          "Any questions please email Tim Nieuwenhuis at:
                   tnieuwe1@jhmi.edu"
        ),
        # Output table
        # mainPanel(dataTableOutput("table"))
        mainPanel(
          tabsetPanel(
            tabPanel("HPAStainR Output", dataTableOutput("table")),
            tabPanel("Summary of Proteins for Cell Types",
                     dataTableOutput("summ_tab"))
          )
        )
      )
    ), # Newly required comma as shiny app is no longer at the end of the data






    server = function(input, output) {



      # Rename cell type data file

      if (is.null(cell_type_data)) {
        cell_type_data <- data.frame(No_data = paste(
                                       "No data given, please use",
                                       "hpa_summary_maker() to make the",
                                       "summary if you wish to see this data",
                                       "in the shiny app")
        )
      } else {
        cell_type_data <- cell_type_data %>% rename(
          `Cell Type` = tissue_cell,
          `Proteins Tested` = proteins,
          `Proteins Detected` = detected,
          `Ratio of Detected to Tested` = det_o_test
        )
      }


      # Insert required function


      # Removed because HPAStainR is loaded in the space?





      output$gene_list <- renderText({
        input$gene_list
      })

      n1 <- reactive({
        (
          HPAStainR(
            gene_list = input$gene_list,
            # gene_list = unlist(str_split(input$gene_list), ','),
            hpa_dat = hpa_dat,
            cancer_dat = cancer_dat,
            tissue_level = ifelse(input$tissue_level == "Yes", TRUE, FALSE),
            stringency = input$stringency,
            scale_abundance = as.logical(input$scale_abundance),
            round_to = input$round_to,
            csv_names = as.logical(input$csv_names),
            percent_or_count = input$percent_or_count,
            stained_gene_data = as.logical(input$stain_gene_results),
            tested_protein_column = as.logical(input$tested_protein_column),
            cancer_analysis = input$cancer_analysis,
            drop_na_row = as.logical(input$drop_na_rows_in),
            adjusted_pvals = as.logical(input$adjusted_pvals)
          )


        )
      })




      output$table <- renderDataTable(print(n1()))


      output$summ_tab <- renderDataTable(print(cell_type_data))


      output$downloadData <- downloadHandler(
        filename = "HPAstainR_results.csv",
        content = function(file) {
          write.csv(print(n1()), file, row.names = FALSE)
        }
      )
    }
  )
}
