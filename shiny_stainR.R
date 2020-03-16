library(shiny)

library(tidyverse)

ui <- fluidPage(
  
  titlePanel("HPAStain.R"),
  #Option selection 
  sidebarLayout(  
    #List of genes
  sidebarPanel(textAreaInput("gene_list", "List of Genes (comma seperated):", "PRSS1,PNLIP,CELA3A,PRL"),
  selectInput("tissue_level", "Cells broken down by tissue site?", c("Yes", "No")),
  selectInput("percent_or_count", "Do you want a count of how often the protein stains in a tissue and/or what percent of your list
              stains in a given cell type and at what level?", c("percent", "count","both")),
               #This scales the genes fo the percentage so you don't get 100 enrichment scores for tissues where 1 protein was tested
  selectInput("scale_abundance", "Scale results for genes that have available data", c("True", "False")),
  selectInput("stringency", "Confidence level of HPA data", c("Normal", "High", "Low"), selected = "Normal"),
  numericInput("round_to", "Round values to:", value =2,  min = 2, max = 5 ),
  selectInput("csv_names", "Column names easier to deal with in csv format", c("False", "True")),
  submitButton(text = "Run HPAStain.R"),
  textOutput("gene_list"),
  downloadButton("downloadData", "Download Output Table")),
  #Output table
  mainPanel(DT::dataTableOutput("table"))
  
  )

)


  #Insert required function
  
  stainR <- function(gene_list, hpa_dat,# weight_stain = F, weight_reliability = F,
                     tissue_level = T,
                     stringency = "normal",
                     scale_abundance= T,
                     scale_genes = T, #Do not include in shiny, this is for my personal pipeline
                     round_to = 2,
                     #subset_genes = T, #unused
                     csv_names= F, 
                     percent_or_count = c("percent", "count", "both")){
   
    #Make gene list robust to incongruencies
    gene_list = gsub(pattern = " ", replacement =  "", x =  gene_list)
    gene_list = unlist(str_split(gene_list, ','))
    gene_list = toupper(gene_list)
    p_o_c = percent_or_count[1]
    
    if (tissue_level == T) {
      cell_o_tiss = "tissue_cell"
      
    }else{
      cell_o_tiss = "Cell.type"
    }
    
    
    
    #Remove any blanks from the gene list
    gene_list <- gene_list[gene_list != ""]
    
    
    #make new column combining cell type and tissue
    hpa_dat <- hpa_dat %>% 
      mutate(Tissue = gsub('[[:digit:]]+', '', Tissue), tissue_cell = paste0(toupper(str_trim(Tissue))," - ", Cell.type))
    
    
    #Set up all cell types for later
    all_cell_types = unique(hpa_dat[[cell_o_tiss]])
    
    #Subset just to genes of interest 
    sub_dat <- subset(hpa_dat, hpa_dat$Gene.name %in% gene_list)
    
    #Remove testis as their cell over stain and 
    sub_dat <- sub_dat %>% filter(Tissue != "testis", !is.na(Cell.type)) 
    
    #Below selects the tolerance of bad data, I suggest normal or high
    if (stringency == "normal") {
      sub_dat <- subset(sub_dat, sub_dat$Reliability %in% c("Enhanced", "Supported")) 
    }
    
    if (stringency == "high") {
      sub_dat <- subset(sub_dat, sub_dat$Reliability %in% c("Enhanced")) 
    }
    
    
    
    
    #Test how many genes are in hpa
    percent_coding <- sum(gene_list %in% sub_dat$Gene.name)/length(gene_list)
    
    #What are the genes in the dataset
    prot_genes <- gene_list[(gene_list %in% sub_dat$Gene.name)]
    
    #What genes are not in the dataset
    non_coding <- gene_list[!(gene_list %in% sub_dat$Gene.name)]
    
    
    #Below code returns dataframe full of NAs in the case of no protein coding genes, this way you know where data is missing
    if (percent_coding == 0) {
      no_dat_matrix  <- matrix(data = NA, nrow = length(all_cell_types), ncol = 7)
      rownames(no_dat_matrix) <- all_cell_types  
      colnames(no_dat_matrix) <- c("High", "Medium", "Low", "Not detected", "enriched_score", "num_genes", "genes")
      
      
      no_dat_tib <-  as_tibble(no_dat_matrix, rownames = "cell_type" )
      
      return(no_dat_tib)
    }
    #Move above section to the end so it can react to all data
    
    
    
    #CELL TYPE ENRICHMENT
    
    #Find levels of expression in cell types
    cell_type_dat  <- table(sub_dat[[cell_o_tiss]], sub_dat$Level)
    
    cell_type_dat_df <- as.data.frame.matrix(cell_type_dat) 
    
    #WE GOTTA scale for tissues known
    #tiss_scale <- rowSums((table(sub_dat$Cell.type, sub_dat$Tissue))/length(prot_genes))
    
    #cell_type_dat_tiss_scale = cell_type_dat/tiss_scale
    
    #Normalize based on how many times they are detected
    cell_type_dat_per <-   apply(cell_type_dat, 2, FUN = function(x){x/rowSums(cell_type_dat)})
    
    
    scaled_for_genes  <- rowSums(table(sub_dat[[cell_o_tiss]],sub_dat$Gene.name))/length(prot_genes)
    
    

    
    
    
    #Scale for only a few genes existing
    if (scale_abundance == T) {
      #Get unique tissue counts
      uni_tiss  <- sub_dat %>% select((!!sym(cell_o_tiss)), Tissue) %>% distinct %>% arrange((!!sym(cell_o_tiss)))
      
      scaled_4_tiss_n_genes  <- scaled_for_genes/ rowSums(table(uni_tiss[[cell_o_tiss]], uni_tiss$Tissue))
      
      
      cell_type_dat_per <- cell_type_dat_per * scaled_4_tiss_n_genes
      
      #New section that fixes counts of differnt groups
      tiss_cell_percent <- sub_dat %>% select(Gene.name,Cell.type, tissue_cell) 
      
      group_list  <- table(tiss_cell_percent[[cell_o_tiss]])
      
      unique_counts <- sapply(group_list, function(x){ifelse(x/length(prot_genes) > 1, (x/length(prot_genes)), 1 )})
      
      
      cell_type_dat_per <- cell_type_dat_per/unique_counts
      
    }
    
    
    
    
    #Below add column if the low medium or high columns don't exist
    if (!("Low" %in% colnames(cell_type_dat_per))) {
      Low <- matrix(0,nrow = nrow(cell_type_dat_per))
      cell_type_dat_per <- cbind(cell_type_dat_per, Low)
      colnames(cell_type_dat_per)[ncol(cell_type_dat_per)] <- "Low"
      #DF
      Low <- matrix(0,nrow = nrow(cell_type_dat_df))
      cell_type_dat_df <- cbind(cell_type_dat_df, Low)
      colnames(cell_type_dat_df)[ncol(cell_type_dat_df)] <- "Low"
      
      
      rm(Low)
    }
    
    if (!("Medium" %in% colnames(cell_type_dat_per))) {
      Medium <- matrix(0,nrow = nrow(cell_type_dat_per))
      cell_type_dat_per <- cbind(cell_type_dat_per, Medium)
      colnames(cell_type_dat_per)[ncol(cell_type_dat_per)] <- "Medium"
      
      #DF
      Medium <- matrix(0,nrow = nrow(cell_type_dat_df))
      cell_type_dat_df <- cbind(cell_type_dat_df, Medium)
      colnames(cell_type_dat_df)[ncol(cell_type_dat_df)] <- "Medium"
      
      
      rm(Medium)
    }
    if (!("High" %in% colnames(cell_type_dat_per))) {
      High <- matrix(0,nrow = nrow(cell_type_dat_per))
      cell_type_dat_per <- cbind(cell_type_dat_per, High)
      colnames(cell_type_dat_per)[ncol(cell_type_dat_per)] <- "High"
      #DF
      High <- matrix(0,nrow = nrow(cell_type_dat_df))
      cell_type_dat_df <- cbind(cell_type_dat_df, High)
      colnames(cell_type_dat_df)[ncol(cell_type_dat_df)] <- "High"
      
      rm(High)
    }

    
        
    #Add all missing cell types, this allows you to see what cells there is no information for
    #This is done below by creating a matrix of NAs of the not included cells
    not_included_cells <- all_cell_types[!(all_cell_types %in% row.names(cell_type_dat_per))]
    
    not_included_matrix  <- matrix(data = NA, nrow = length(not_included_cells), ncol = ncol(cell_type_dat_per))
    
    rownames(not_included_matrix) <- not_included_cells
    
    cell_type_dat_per <- rbind(cell_type_dat_per, not_included_matrix)
    
    cell_type_dat_mat <- rbind(as.matrix(cell_type_dat_df), not_included_matrix)
    
    
    
    #Below is where we generate the final tibble
    cell_type_out <-  as_tibble(cell_type_dat_per, rownames = "cell_type" ) %>% #1. Make the data a tibble with rownames as cell_type
      #dplyr::arrange(desc(High), desc(Medium), desc(Low), desc(`Not detected`)) %>% #2. No longer used step
      mutate_if(is.numeric, round, round_to) %>% #3. 
      dplyr::select(cell_type, High, Medium, Low, `Not detected`) %>% #4
      mutate(enriched_score = (High * 100) + (Medium * 50) + (Low *25),
             num_genes =  sum(gene_list %in% sub_dat$Gene.name)) %>%
      arrange(desc(enriched_score)) 
    
    
    
    #Prepare count data for joining
    cell_type_count <- as_tibble(cell_type_dat_mat, rownames = "cell_type") %>% rename("high_expression_count" = "High", 
                                                                                       "medium_expression_count" = "Medium",
                                                                                       "low_expression_count" = "Low",
                                                                                       "not_detcted_count" = "Not detected")
    
    
    cell_type_out <- left_join(cell_type_out, cell_type_count)
    
    
   
    
    #Change genes in column to only thos detcted
    if (scale_genes == T) {
      prot_genes
      tiss_gene_table  <- table(sub_dat[[cell_o_tiss]],sub_dat$Gene.name) > 0.5
      cell_types_current <- cell_type_out$cell_type
      cell_types_current
      
      
      
      
      
      gene_col  <- NULL
      gene_count <- NULL
      for (cells_in in cell_types_current) {
        
        #Add grouped or split gene option here
        
        
        #cells_in <- "chondrocytes"
        if (cells_in %in% rownames(tiss_gene_table)) {
          temp_genes <- names(tiss_gene_table[cells_in,])[tiss_gene_table[cells_in,] == T]
          gene_col <-  c(gene_col, paste0(temp_genes,collapse = ", "))
          gene_count <- c(gene_count, length(temp_genes))
          
        }else{
          gene_col <- c(gene_col, "")
          gene_count <- c(gene_count, 0)
        }
        
        
        
      }
      
    }
    
    if (scale_genes == T) {
      
      cell_type_out$genes <- as.vector(gene_col)
      cell_type_out$num_genes <- gene_count
    }else{
      cell_type_out$genes <- paste0(prot_genes,collapse = ", ")
    }
    
   
    #Change names; might need to change once count data is incorporated
    if (csv_names == T) {
      cell_type_out <- cell_type_out %>% select(cell_type,
                                                percent_high_expression = High,
                                                high_expression_count,
                                                percent_medium_expression = Medium,
                                                medium_expression_count,
                                                percent_low_expression = Low,
                                                low_expression_count,
                                                percent_not_detected = `Not detected`,
                                                not_detcted_count,
                                                number_of_genes = num_genes,
                                                everything())
    }
    
    if (csv_names == F) {
      cell_type_out <- cell_type_out %>% select(`Cell Type` = cell_type,
                                                `Percent High Expression` = High,
                                                `High Expression Count` = high_expression_count,
                                                `Percent Medium Expression` = Medium, 
                                                `Medium Expression Count` = medium_expression_count,
                                                `Percent Low Expression`= Low,
                                                `low_expression_count` = low_expression_count,
                                                `Percent Not Detected` = `Not detected`,
                                                `Not Detected Count` = not_detcted_count, 
                                                `Number of Genes` = num_genes,
                                                `Enriched Score` = enriched_score,
                                                Genes = genes,
                                                everything())
    }
    
    #Select count percent prefernce
    if (p_o_c == "percent") {
      cell_type_out <- cell_type_out[, -grep("ount", colnames(cell_type_out))]
      
    }
    
    if (p_o_c == "count") {
      cell_type_out <- cell_type_out[, -grep("ercent", colnames(cell_type_out))]
      
    }
    
    
    return((cell_type_out))
    
  }
  
  
  
  


server <- function(input, output){
  output$gene_list <- renderText({input$gene_list})
  
  n1 <- reactive({
    (stainR(
      gene_list =input$gene_list,
      #gene_list = unlist(str_split(input$gene_list), ','),
      hpa_dat = hpa_dat, 
      tissue_level = ifelse(input$tissue_level == "Yes", T, F),
      stringency = input$stringency,
      scale_abundance = as.logical(input$scale_abundance),
      round_to = input$round_to,
      csv_names = as.logical(input$csv_names),
      percent_or_count = input$percent_or_count
      
      
      
    ))
    
    
  })
  
  output$table <- DT::renderDataTable(print(n1()))

    
  
 
 
  output$downloadData <- downloadHandler(
    filename ="HPAstainR_results.csv",
    content = function(file){
      write.csv(print(n1()), file, row.names = F)
    }
      )
  

  
  
  
  
}

shinyApp(ui, server)

