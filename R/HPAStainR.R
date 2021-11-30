#' @title HPAStainR
#'
#' @description Uses a protein/gene list to query Human Protein Atlas (HPA)
#'    staining data.
#'
#' @param gene_list A list of proteins or genes that you want to query the HPA
#'    staining data with.
#' @param hpa_dat The data frame of normal HPA staining data data, required to
#'    run HPAStainR.
#' @param cancer_dat The data frame of pathologic HPA staining data, required
#'    to run HPAStainr.
#' @param cancer_analysis A character string indicating inclusion of cancer data
#'    in the result, must be one of 'normal' (default), 'cancer', or 'both'.
#' @param tissue_level A boolean that determines whether tissue level data for
#'    the cell types are included. Default is TRUE
#' @param stringency A character string indicating how stringent the confidence
#'    level of the staining findings have to be. Must be 'normal' (default),
#'    'high', or 'low'. This stringency is based on the `Reliability` column
#'    from the hpa_dat object which varies from "Enhanced", "Supported",
#'    "Approved", to "Uncertain" in decreasing order of certainty. Low
#'    stringency includes all data, normal stringency includes "Enhanced",
#'    "Supported", and "Approved", while high stringency only includes
#'    "Enhanced" and "Supported". Further information about these
#'     categorizations can be found in the following link
#'     https://www.proteinatlas.org/about/assays+annotation
#' @param scale_abundance A boolean that determines whether you scale Staining
#'    Score based on the size of the gene list. Default is TRUE.
#' @param round_to  A numeric that determines how many decimals in numeric
#'    outputs are desired. Default 2.
#' @param csv_names A Boolean determining if you want names suited for a csv
#'    file/pipeline, or for presentation. Default is TRUE giving csv names.
#' @param stained_gene_data A boolean determining if there is a list of which
#'    proteins stained, TRUE is default.
#' @param tested_protein_column A boolean determining if there is a column
#'    listing which proteins were tested, TRUE is default.
#' @param percent_or_count A character string determining if percent of proteins
#'    stained, count of proteins stained, or both are shown for high, medium,
#'     and low staining. Must be 'percent' (default), 'count', or 'both'.
#' @param drop_na_row A boolean that determines if cell types with no proteins
#'    tested are kept or dropped, default is FALSE.
#' @param test_type A character vector for either "fisher" or "chi square",
#'    used to select the statistical test for determining cell type enrichment.
#'    The two options are Fisher's Exact Test and a Chi Square test. The
#'    original version of HPAStainR defaulted to the Chi Square test, however
#'    because this requires simulated values to run correctly, we suggest the
#'    usage of the Fisher's Exact Test for consistency.
#' @param adjusted_pvals A boolean indicating if you want the p-values corrected
#'    for multiple testing. Default is TRUE.
#'
#' @section Details:
#' Calculation of the staining score below:
#' \deqn{(\frac{h \times 100}{t}) + (\frac{m \times 50}{t}) +
#'  (\frac{l \times 25}{t})}
#'
#'
#' @return    A tibble containing the results of HPAStainR.
#'
#' @examples
#'    ## Below will give you the results found on the shiny app website
#'    ## This example also uses HPA_data_downloader output as an example
#'    HPA_data <- HPA_data_downloader(tissue_type = 'both', save_file = FALSE)
#'    HPA_out <- HPAStainR(c('PRSS1', 'PNLIP', 'CELA3A', 'PRL'),
#'    HPA_data$hpa_dat,
#'    HPA_data$cancer_dat,
#'    'both')
#' @import dplyr 
#' @import shiny
#' @import tibble
#' @import tidyr
#' @importFrom scales percent
#' @importFrom stringr str_detect str_split str_trim
#' @importFrom stats chisq.test p.adjust quantile fisher.test
#' @export

HPAStainR <- function(gene_list,
                      hpa_dat, cancer_dat = data.frame(),
                      cancer_analysis = c("normal", "cancer", "both"),
                      tissue_level = TRUE,
                      stringency = c("normal", "high", "low"),
                      scale_abundance = TRUE,
                      round_to = 2,
                      csv_names = TRUE,
                      stained_gene_data = TRUE,
                      tested_protein_column = TRUE,
                      percent_or_count = c("percent", "count", "both"),
                      drop_na_row = FALSE,
                      test_type = c("fisher", "chi square"),
                      adjusted_pvals = TRUE) {
    ### Hi Maisa, so this is the main function of HPAStainR, HPAStainR! Our
    ### overall goal for this "rotation" you're doing here is to make the code
    ### faster, while maintaining the current outputs. My advice would be to
    ### see what areas of the code could be placed into nested functions
    ### (turning these blocks of code into functions will make future changes,
    ###  unit testing, and speed testing easier as well). I will do my best to
    ### assist you in this project, but obviously this will also required your
    ### own initiative and problem solving. However if you help me here I will
    ### make you a co-author on this package will will be great for your future
    ### career in my opinion! Thanks for doing this and always feel free to ping
    ### me!
    
    ### Catch issues
    ## Is the HPA or Cancer dataframe ther or working
    cancer_analysis <- cancer_analysis[1]
    if (cancer_analysis == "normal" | cancer_analysis == "both") {
        if (!is.data.frame(hpa_dat)) {
            stop(paste0("Are you sure hpa_dat is a dataframe? Download through",
                       "HPADownloader if you're having issues or use hpar"))
        }
    }
    
    if (cancer_analysis == "cancer" | cancer_analysis == "both") {
        if (!is.data.frame(cancer_dat)) {
            stop(paste("Are you sure cancer_dat is a dataframe? Download",
                       "through HPADownloader if you're having issues"))
        }
    }
    ## In case no one select an option this will pick the default
    stringency <- stringency[1]
    ### Maisa Note:
    ## This was aneasy way to make cancer only analysis work, though inefficient, as it results
    ## in both the cancer analysis and normal analysis running in tandem
    ## this could be replaced probably through using a list and lapply.
    cancer_only <- FALSE
    if (cancer_analysis == "cancer") {
        cancer_analysis <- "both"
        cancer_only <- TRUE
    }
    ### Maisa Note:
    ## A holdover from my personal pipeline
    ## this should be removed if possible
    scale_genes <- TRUE
    
    ## Select the correct test type
    test_type = test_type[1]
    
    ### Make gene list robust to different inputs test if comma 
    ## separated or non comma separated
    
    if (!(str_detect(gene_list, ","))[1]) {
        gene_list <- gsub(pattern = "\\s", replacement = ",", x = gene_list)
        gene_list <- gsub(pattern = ",{2,}", replacement = ",", x = gene_list)
    }
    
    gene_list <- gsub(pattern = " ", replacement = "", x = gene_list)
    gene_list <- unlist(str_split(gene_list, ","))
    gene_list <- toupper(gene_list)
    p_o_c <- percent_or_count[1]
    
    ## This feels redundant as tissue level could replace the cell_o_tiss
    ## object, but removing it is not high priority.
    if (tissue_level == TRUE) {
        cell_o_tiss <- "tissue_cell"
    } else {
        cell_o_tiss <- "Cell.type"
    }
    
    
    ## Remove any blanks from the gene list
    gene_list <- gene_list[gene_list != ""]
    
    
    ## Make new column combining cell type and tissue
    hpa_dat <- hpa_dat %>%
        mutate(Tissue = gsub("[[:digit:]]+", "", Tissue),
               tissue_cell = paste0(toupper(str_trim(Tissue)),
                                    " - ",
                                    Cell.type))
    
    
    ## Set up all cell types for later
    all_cell_types <- unique(hpa_dat[[cell_o_tiss]])
    
    ## Subset just to genes of interest
    sub_dat <- subset(hpa_dat, hpa_dat$Gene.name %in% gene_list)
    
    # Remove testis as their cells over stain and remove NA cells
    sub_dat <- sub_dat %>% filter(Tissue != "testis", !is.na(Cell.type),
                                  tissue_cell != "N/A - N/A")
    
    # Below selects the tolerance of bad data, I suggest normal or high
    if (stringency == "normal") {
        sub_dat <- subset(sub_dat, sub_dat$Reliability %in%
                              c("Enhanced", "Supported", "Approved"))
    }
    
    if (stringency == "high") {
        sub_dat <- subset(sub_dat, sub_dat$Reliability %in%
                              c("Enhanced", "Supported"))
    }
    
    ## Test how many genes are in hpa
    percent_coding <- sum(gene_list %in% sub_dat$Gene.name)/length(gene_list)
    ## What are the genes in the dataset
    prot_genes <- gene_list[(gene_list %in% sub_dat$Gene.name)]
    ## What genes are not in the dataset
    non_coding <- gene_list[!(gene_list %in% sub_dat$Gene.name)]
    ## Below code returns dataframe full of NAs in the case of no protein 
    ## coding genes, this way you know where data is missing
    if (percent_coding == 0) {
        no_dat_matrix <- matrix(data = NA, nrow = length(all_cell_types),
                                ncol = 7)
        rownames(no_dat_matrix) <- all_cell_types
        colnames(no_dat_matrix) <- c("High",
                                     "Medium",
                                     "Low",
                                     "Not detected",
                                     "staining_score",
                                     "num_genes",
                                     "genes")
        
        no_dat_tib <- as_tibble(no_dat_matrix, rownames = "cell_type")
        return(no_dat_tib)
    }
    ## CELL TYPE ENRICHMENT------------ Find levels of expression in cell types
    cell_type_dat <- table(sub_dat[[cell_o_tiss]], sub_dat$Level)
    cell_type_dat_df <- as.data.frame.matrix(cell_type_dat)
    ## cell_type_dat_tiss_scale = cell_type_dat/tiss_scale Normalize based on
    ## how many times they are detected
    cell_type_dat_per <- apply(cell_type_dat, 2, FUN = function(x) {
        x/rowSums(cell_type_dat)
    })
    scaled_for_genes <- rowSums(table(sub_dat[[cell_o_tiss]],
                                      sub_dat$Gene.name))/length(prot_genes)
    ## Scale for only a few genes existing
    if (scale_abundance == TRUE) {
        ## Get unique tissue counts
        uni_tiss <- sub_dat %>%
            select((!!sym(cell_o_tiss)), Tissue) %>%
            distinct() %>% arrange((!!sym(cell_o_tiss)))
        scaled_4_tiss_n_genes <- scaled_for_genes/
            rowSums(table(uni_tiss[[cell_o_tiss]], uni_tiss$Tissue))
        cell_type_dat_per <- cell_type_dat_per * scaled_4_tiss_n_genes
        ## New section that fixes counts of different groups
        tiss_cell_percent <- sub_dat %>% select(Gene.name, Cell.type,
                                                tissue_cell)
        group_list <- table(tiss_cell_percent[[cell_o_tiss]])
        unique_counts <- vapply(group_list, function(x) {
            ifelse(x/length(prot_genes) > 1, (x/length(prot_genes)), 1)
        }, numeric(1))
        cell_type_dat_per <- cell_type_dat_per/unique_counts
    }
    
    ### cell_type_dat_per unit test here
    ### cell_type_dat_per_A
    
    ## Below add column if the low medium or high columns don't exist--------
    ### Maisa Note: Obviously this copy and paste malarky below is bad code
    ### we should try to tighten it up.
    if (!("Low" %in% colnames(cell_type_dat_per))) {
        Low <- matrix(0, nrow = nrow(cell_type_dat_per))
        cell_type_dat_per <- cbind(cell_type_dat_per, Low)
        colnames(cell_type_dat_per)[ncol(cell_type_dat_per)] <- "Low"
        ## DF
        Low <- matrix(0, nrow = nrow(cell_type_dat_df))
        cell_type_dat_df <- cbind(cell_type_dat_df, Low)
        colnames(cell_type_dat_df)[ncol(cell_type_dat_df)] <- "Low"
        
        rm(Low)
    }
    
    if (!("Medium" %in% colnames(cell_type_dat_per))) {
        Medium <- matrix(0, nrow = nrow(cell_type_dat_per))
        cell_type_dat_per <- cbind(cell_type_dat_per, Medium)
        colnames(cell_type_dat_per)[ncol(cell_type_dat_per)] <- "Medium"
        
        ## DF
        Medium <- matrix(0, nrow = nrow(cell_type_dat_df))
        cell_type_dat_df <- cbind(cell_type_dat_df, Medium)
        colnames(cell_type_dat_df)[ncol(cell_type_dat_df)] <- "Medium"
        
        
        rm(Medium)
    }
    if (!("High" %in% colnames(cell_type_dat_per))) {
        High <- matrix(0, nrow = nrow(cell_type_dat_per))
        cell_type_dat_per <- cbind(cell_type_dat_per, High)
        colnames(cell_type_dat_per)[ncol(cell_type_dat_per)] <- "High"
        ## DF
        High <- matrix(0, nrow = nrow(cell_type_dat_df))
        cell_type_dat_df <- cbind(cell_type_dat_df, High)
        colnames(cell_type_dat_df)[ncol(cell_type_dat_df)] <- "High"
        
        rm(High)
    }
    
    ### cell_type_dat_per unit test here
    ### cell_type_dat_per_B
    
    ## Add all missing cell types, this allows you to see what cells there is no
    ## information for--------
    
    ## This is done below by creating a matrix of NAs of the not included cells
    not_included_cells <- all_cell_types[!(all_cell_types %in%
                                               row.names(cell_type_dat_per))]
    
    not_included_matrix <- matrix(data = NA,
                                  nrow = length(not_included_cells),
                                  ncol = ncol(cell_type_dat_per))
    
    rownames(not_included_matrix) <- not_included_cells
    
    cell_type_dat_per <- rbind(cell_type_dat_per, not_included_matrix)
    
    cell_type_dat_mat <- rbind(as.matrix(cell_type_dat_df), not_included_matrix)
    
    ### Maisa Note: Overall Cancer code needs to be cleaned up some way so that
    ### it isn't completely separate from the normal tissue code
    ## Adding CANCER dat
    if (cancer_analysis == "both" | cancer_analysis == "Only") {
        sub_cancer <- cancer_dat %>% filter(Gene.name %in% gene_list)
        
        sub_cancer <- sub_cancer %>% filter(!is.na(Cancer), !is.na(High))
        
        cancer_count <- sub_cancer %>%
            group_by(Cancer) %>%
            summarise(High = sum(High),
                      Medium = sum(Medium),
                      Low = sum(Low),
                      `Not detected` = sum(Not.detected))
        ## Remove rownames due to dplyr
        rownames(cancer_count) <- c()
        
        cancer_per <- cancer_count
        
        cancer_per[, -1] <- (cancer_count[, -1]/rowSums(cancer_count[, -1]))
        
        cancer_count <- as.matrix(cancer_count %>%
                                      column_to_rownames(var = "Cancer"))
        
        cancer_per <- as.matrix(cancer_per %>%
                                    column_to_rownames(var = "Cancer"))
        
        
        if (cancer_analysis == "both") {
            cell_type_dat_mat <- rbind(cell_type_dat_mat, cancer_count)
            cell_type_dat_per <- rbind(cell_type_dat_per, cancer_per)
        }
    }
    
    
    ### cell_type_dat_per unit test here
    ### cell_type_dat_per_C
    
    ## Before tibble creation catch and fix lack of `Not detected` column
    if (!("Not detected" %in% colnames(cell_type_dat_per))) {
        na_ind <- is.na(apply(cell_type_dat_per,1,sum))
        `Not detected` <- ifelse(na_ind == TRUE, NA, 0)
        cell_type_dat_per <- cbind(cell_type_dat_per, `Not detected`)
        cell_type_dat_mat <- cbind(cell_type_dat_mat, `Not detected`)
    }
    
    
    
    ## Below is where we generate the final tibble----------------- 
    ## 1. Make the data a tibble with rownames as
    ## cell_type
    cell_type_out <- as_tibble(cell_type_dat_per, rownames = "cell_type") %>%
        # 2. Round all numeric columns according to arguement
    mutate_if(is.numeric, round, round_to) %>%
        # 3. Select columns
    select(cell_type, High, Medium, Low, `Not detected`) %>%
        # 4. Generate the staining score
    mutate(staining_score = (High * 100) + (Medium * 50) + (Low * 25),
           num_genes = sum(gene_list %in% sub_dat$Gene.name)) %>% 
        arrange(desc(staining_score))

    
    ## Prepare count data for joining
    cell_type_count <- as_tibble(cell_type_dat_mat, rownames = "cell_type") %>%
        rename(high_expression_count = "High",
               medium_expression_count = "Medium",
               low_expression_count = "Low",
               not_detected_count = "Not detected")
    
    
    cell_type_out <- left_join(cell_type_out,
                               cell_type_count,
                               by = "cell_type")
    
    ### cell_type_out unit test here
    ### cell_type_out_A    
    
    
    ## Change genes in column to only those detected-------------
    
    
    if (scale_genes == TRUE) {
        prot_genes
        tiss_gene_table <- table(sub_dat[[cell_o_tiss]],
                                 sub_dat$Gene.name) > 0.5
        ## CANCER
        if (cancer_analysis == "both") {
            cancer_gene_table <- as.matrix(table(sub_cancer$Cancer,
                                                 sub_cancer$Gene.name) > 0.5)
            ## Make it robust for non matching
            if (ncol(tiss_gene_table) != ncol(cancer_gene_table)) {
                not_shared_normal <- colnames(
                    cancer_gene_table)[!(colnames(cancer_gene_table) %in%
                                             colnames(tiss_gene_table))]
                not_shared_cancer <- colnames(
                    tiss_gene_table)[!(colnames(tiss_gene_table) %in%
                                           colnames(cancer_gene_table))]
                
                norm_add_matrix <- matrix(
                    data = FALSE, nrow = nrow(tiss_gene_table),
                    ncol = length(not_shared_normal))
                colnames(norm_add_matrix) <- not_shared_normal
                tiss_gene_table <- (cbind(tiss_gene_table, norm_add_matrix))
                
                
                cancer_add_matrix <- matrix(
                    data = FALSE, nrow = nrow(cancer_gene_table),
                    ncol = length(not_shared_cancer))
                colnames(cancer_add_matrix) <- not_shared_cancer
                cancer_gene_table <- (cbind(cancer_gene_table,
                                            cancer_add_matrix))
            }
            
            tiss_gene_table <- rbind(tiss_gene_table, cancer_gene_table)
        }
        
        cell_types_current <- cell_type_out$cell_type
        
        gene_col <- NULL
        gene_count <- NULL
        for (cells_in in cell_types_current) {
            
            ## Add grouped or split gene option here
            if (cells_in %in% rownames(tiss_gene_table)) {
                temp_genes <- names(
                    tiss_gene_table[cells_in, ])[tiss_gene_table[cells_in,
                                                                 ] == TRUE]
                gene_col <- c(gene_col, paste0(temp_genes, collapse = ", "))
                gene_count <- c(gene_count, length(temp_genes))
            } else {
                gene_col <- c(gene_col, "")
                gene_count <- c(gene_count, 0)
            }
        }
    }
    
    if (scale_genes == TRUE) {
        cell_type_out$genes <- as.vector(gene_col)
        cell_type_out$num_genes <- gene_count
    } else {
        cell_type_out$genes <- paste0(prot_genes, collapse = ", ")
    }
    
    
    ## Add the option that gives a column if a gene is available
    
    ## Tissue level determines if the analysis looks at cell types on a tissue
    ## level instead of a general
    ## level------
    
    ## This is ideal as cell types in different tissues, though similar, have
    ## different profiles.
    if (tissue_level == TRUE) {
        staining_dat <- sub_dat %>% filter(Level != "Not detected")
        staining_tf_df <- as.matrix.data.frame(
            table(staining_dat$tissue_cell, staining_dat$Gene.name) > 0, TRUE)
        colnames(staining_tf_df) <- colnames(
            table(staining_dat$tissue_cell, staining_dat$Gene.name))
        staining_tf_df <- as.data.frame(staining_tf_df)
        
        ## CANCER
        if (cancer_analysis == "both" | cancer_analysis == "Only") {
            cancer_staining_dat <- sub_cancer %>%
                filter(High > 0 | Low > 0 | Medium > 0)
            cancer_staining_tf_df <- as.matrix.data.frame(
                table(cancer_staining_dat$Cancer,
                      cancer_staining_dat$Gene.name) > 0, TRUE)
            colnames(cancer_staining_tf_df) <- colnames(
                table(cancer_staining_dat$Cancer,
                      cancer_staining_dat$Gene.name))
            cancer_staining_tf_df <- as.data.frame(cancer_staining_tf_df)
            ## Put in blank information
            false_cols <- colnames(
                staining_tf_df[!(colnames(staining_tf_df) %in%
                                     colnames(cancer_staining_tf_df))])
            ## Make false_matrix
            false_matrix <- matrix(data = FALSE,
                                   ncol = length(false_cols),
                                   nrow = nrow(cancer_staining_tf_df))
            colnames(false_matrix) <- false_cols
            cancer_staining_tf_df <- cbind(cancer_staining_tf_df, false_matrix)
            which(colnames(cancer_staining_tf_df) %in% colnames(staining_tf_df))
            ## Reorder
            new_order <- match(colnames(staining_tf_df),
                               colnames(cancer_staining_tf_df))
            cancer_staining_tf_df <- cancer_staining_tf_df[, new_order]
            
            if (cancer_analysis == "both") {
                staining_tf_df <- rbind(staining_tf_df, cancer_staining_tf_df)
            }
        }
        ## Cancer end
    } else {
        staining_dat <- sub_dat %>% filter(Level != "Not detected")
        staining_tf_df <- as.matrix.data.frame(
            table(staining_dat$Cell.type,
                  staining_dat$Gene.name) > 0, TRUE)
        colnames(staining_tf_df) <- colnames(
            table(staining_dat$Cell.type, staining_dat$Gene.name))
        staining_tf_df <- as.data.frame(staining_tf_df)
        
        ## CANCER
        if (cancer_analysis == "both" | cancer_analysis == "Only") {
            cancer_staining_dat <- sub_cancer %>% filter(High > 0 | Low > 0 |
                                                             Medium > 0)
            cancer_staining_tf_df <- as.matrix.data.frame(
                table(cancer_staining_dat$Cancer,
                      cancer_staining_dat$Gene.name) > 0, TRUE)
            colnames(cancer_staining_tf_df) <- colnames(
                table(cancer_staining_dat$Cancer,
                      cancer_staining_dat$Gene.name))
            cancer_staining_tf_df <- as.data.frame(cancer_staining_tf_df)
            ## Put in blank information
            false_cols <- colnames(
                staining_tf_df[!(colnames(staining_tf_df) %in%
                                     colnames(cancer_staining_tf_df))])
            ## Make false_matrix
            false_matrix <- matrix(data = FALSE,
                                   ncol = length(false_cols),
                                   nrow = nrow(cancer_staining_tf_df))
            colnames(false_matrix) <- false_cols
            cancer_staining_tf_df <- cbind(cancer_staining_tf_df, false_matrix)
            which(colnames(cancer_staining_tf_df) %in% colnames(staining_tf_df))
            ## Reorder
            new_order <- match(colnames(staining_tf_df),
                               colnames(cancer_staining_tf_df))
            cancer_staining_tf_df <- cancer_staining_tf_df[, new_order]
            
            if (cancer_analysis == "both") {
                staining_tf_df <- rbind(staining_tf_df, cancer_staining_tf_df)
            }
        }
        ## cancer end
    }
    
    ## For loop to replace T F with name
    for (col_n in seq_len(ncol(staining_tf_df))) {
        gene <- colnames(staining_tf_df)[col_n]
        staining_tf_df[, col_n] <- ifelse(staining_tf_df[, col_n] == TRUE,
                                          gene, "")
    }
    
    
    stained_list <- apply(as.matrix(staining_tf_df), 1, paste, collapse = ",")
    ## Remove all , at the end of strings
    while (sum(str_detect(string = stained_list, pattern = ",$")) > 0) {
        stained_list <- gsub(pattern = ",$", x = stained_list, replacement = "")
    }
    
    ## Remove all , at beginning
    while (sum(str_detect(string = stained_list, pattern = "^,")) > 0) {
        stained_list <- gsub(pattern = "^,", x = stained_list, replacement = "")
    }
    
    ## Remove all ,, and replace with ,
    while (sum(str_detect(string = stained_list, pattern = ",,")) > 0) {
        stained_list <- gsub(pattern = ",,", x = stained_list,
                             replacement = ",")
    }
    
    stained_list <- gsub(",", ", ", stained_list)
    
    stained_out <- as.data.frame(stained_list, stringsAsFactors = FALSE) %>%
        rownames_to_column(var = "cell_type")
    
    cell_type_out <- left_join(cell_type_out, stained_out, by = "cell_type")
    # }
    
    
    ## Move stained out into upper list
    
    
    
    ## Only cancer
    
    if (cancer_only == TRUE) {
        cell_type_out <- cell_type_out %>% filter(cell_type %in%
                                                      cancer_dat$Cancer)
    }
    
    
    ### CHI TEST HERE ------------- Insert chi square test here, first normal
    ### data, then cancer, then both Make
    ### sure to csv names filter cell type out to remove NAs
    
    
    ## Remove NAs and Testis as we don't use it
    ubi_test <- hpa_dat %>%
        mutate(stained = ifelse(Level != "Not detected", TRUE, FALSE),
               in_list = Gene.name %in% gene_list) %>%
        filter(Tissue != "testis",
               !is.na(Cell.type),
               tissue_cell != "N/A - N/A")
    ## Get a list of tested genes
    genes_tested <- table(ubi_test$Gene.name)
    ## Filter down to stained genes
    ubi_test_filt <- ubi_test %>% filter(stained == TRUE)
    ## Make a table of stained proteins by cell -tisssue
    ubi_table <- table(ubi_test_filt$tissue_cell, ubi_test_filt$Gene.name)
    ## Now just remove non-matching proteins
    ubi_table_filt <- ubi_table[, colnames(ubi_table) %in%
                                    rownames(genes_tested)]
    genes_tested_filt <- genes_tested[rownames(genes_tested) %in%
                                          colnames(ubi_table_filt)]
    ## Create the out object which is a ratio of those stained over those tested
    out <- (colSums(ubi_table_filt)/as.vector(genes_tested_filt))
    ## Cut the data down to rare genes found in less thant the 1st quartile
    ## quart_ind <- out[out <= quantile(out,
    ## .25)]
    quart_ind <- out[out <= quantile(out, 0.15)]
    
    ### The Chi Square calculation Make quart hpa from sub hpa, necessary?
    
    ## testing switching ubi test for sub hpa
    
    
    
    quart_hpa <- ubi_test %>% filter(Gene.name %in% names(quart_ind))
    
    ## Not 2 levels is used to remove cell types that fail to have two levels in
    ## the gene list
    not_2_levels <- quart_hpa %>% group_by(tissue_cell) %>%
        summarise(stain_mean = mean(stained), list_mean = mean(in_list)) %>% 
        filter(stain_mean == 1 |
                   stain_mean == 0 |
                   list_mean == 1 |
                   list_mean == 0)
    
    ## filter down to top specificity genes
    quart_hpa <- quart_hpa %>% filter(!(tissue_cell %in%
                                            not_2_levels$tissue_cell))
    ## The chi test
    
    if (nrow(quart_hpa) != 0) {
        test_out <- quart_hpa %>% group_by(tissue_cell) %>%
            {if (test_type == "chi square") summarise(., p_val = chisq.test(
                stained,
                in_list,
                simulate.p.value = TRUE)$p.value)
                else . } %>% 
            {if (test_type == "fisher") summarise(., p_val = fisher.test(
                stained,
                in_list)$p.value)
                else . } %>%
            rename(cell_type = tissue_cell)
        test_out$p_val_adj <- p.adjust(test_out$p_val)
    }
    
    if (exists("test_out")) {
        if (cancer_analysis == "normal") {
            cell_type_out <- left_join(cell_type_out, test_out, by = "cell_type")
        } else {
            test_out_tiss <- test_out
        }
    }
    
    
    ##### DUPLICATION occurs below
    #### CHI SQUARE CANCER ANALYSIS---------
    if (cancer_analysis == "both" | cancer_analysis == "Only") {
        sub_cell_type <- cell_type_out %>%
            filter(!(is.na(high_expression_count)))
        # 
        sub_canc <- cancer_dat %>% filter(Cancer %in%
                                              (sub_cell_type$cell_type)) %>%
            unique()
        
        sub_canc <- sub_canc %>%
            mutate(stained = ifelse(Not.detected != 0, TRUE, FALSE),
                   in_list = ifelse(Gene.name %in% gene_list, TRUE, FALSE))
        
        ## Remove NAs and Testis as we don't use it
        ubi_test <- cancer_dat %>%
            mutate(stained = ifelse(Not.detected != 0, TRUE, FALSE),
                   in_list = ifelse(Gene.name %in% gene_list, TRUE, FALSE))
        ## Get a list of tested genes
        genes_tested <- table(ubi_test$Gene.name)
        ## Filter down to stained genes
        ubi_test_filt <- ubi_test %>% filter(stained == TRUE)
        ## Make a table of stained proteins by cell -tisssue
        ubi_table <- table(ubi_test_filt$Cancer, ubi_test_filt$Gene.name)
        ## Now just remove non-matching proteins
        ubi_table_filt <- ubi_table[, colnames(ubi_table) %in%
                                        rownames(genes_tested)]
        genes_tested_filt <- genes_tested[rownames(genes_tested) %in%
                                              colnames(ubi_table_filt)]
        ## Create the out object which is a ratio of those stained over those
        ## tested
        out <- (colSums(ubi_table_filt)/as.vector(genes_tested_filt))
        ## Cut the data down to rare genes found in less than the 1st quartile
        quart_ind <- out[out <= quantile(out, 0.15)]
        
        
        quart_canc <- ubi_test %>% filter(Gene.name %in% names(quart_ind))
        
        ## Not 2 levels is used to remove cell types that fail to have two
        ## levels in the gene list
        not_2_levels <- quart_canc %>%
            filter(!is.na(stained)) %>%
            group_by(Cancer) %>%
            summarise(stain_mean = mean(stained), 
                      list_mean = mean(in_list)) %>%
            filter(stain_mean == 1 |
                       stain_mean == 0 |
                       list_mean == 1 | list_mean == 0)
        
        ## filter down to top specificity genes
        quart_canc <- quart_canc %>% filter(!(Cancer %in% not_2_levels$Cancer))
        
        if (nrow(quart_canc) != 0) {
            test_out <- quart_canc %>% group_by(Cancer) %>%
                {if (test_type == "chi square") summarise(., p_val = chisq.test(
                    stained,
                    in_list,
                    simulate.p.value = TRUE)$p.value)
                    else . } %>% 
                {if (test_type == "fisher") summarise(., p_val = fisher.test(
                    stained,
                    in_list)$p.value)
                    else . } %>%
                rename(cell_type = Cancer)
            
            test_out$p_val_adj <- p.adjust(test_out$p_val)
        }
        
        
        
        if (cancer_analysis == "both") {
            test_out <- bind_rows(test_out_tiss, test_out)
            ## For some reason this left_join duplicates rows, and I have no
            ## idea why, but unique fixes it
            cell_type_out <- left_join(cell_type_out,
                                       test_out, by = "cell_type") %>% unique()
        } else {
            cell_type_out <- left_join(cell_type_out,
                                       test_out, by = "cell_type") %>% unique()
        }
    }
    
    
    
    ## Fixing Pvalues------------- In case no pvals
    if (!("p_val" %in% colnames(cell_type_out))) {
        cell_type_out <- cell_type_out %>% mutate(p_val = 1, p_val_adj = 1)
    }
    
    
    ## Change names-----------------
    if (csv_names == TRUE) {
        cell_type_out <- cell_type_out %>%
            select(cell_type,
                    percent_high_expression = High,
                    high_expression_count, 
                    percent_medium_expression = Medium,
                    medium_expression_count,
                    percent_low_expression = Low,
                    low_expression_count, 
                    percent_not_detected = `Not detected`,
                    not_detected_count,
                    number_of_proteins = num_genes,
                    tested_proteins = genes, 
                    detected_proteins = stained_list, everything())
    }
    
    if (csv_names == FALSE) {
        cell_type_out <- cell_type_out %>%
            mutate(High = percent(High, acccuracy = 0.1),
                   Medium = percent(Medium, accuracy = 0.1),
                   Low = percent(Low, accuracy = 0.1),
                   `Not detected` = percent(`Not detected`, accuracy = 0.1)) %>% 
            select(`Cell Type` = cell_type,
                   `Percent High Expression` = High, 
                   `High Expression Count` = high_expression_count,
                   `Percent Medium Expression` = Medium,
                   `Medium Expression Count` = medium_expression_count,
                   `Percent Low Expression` = Low,
                   `Low Expression Count` = low_expression_count,
                   `Percent Not Detected` = `Not detected`,
                   `Not Detected Count` = not_detected_count,
                   `Number of Proteins` = num_genes,
                   `Staining Score` = staining_score,
                   `Tested Proteins` = genes,
                   `Detected Proteins` = stained_list,
                   `P-Value` = p_val,
                   `P-Value Adjusted` = p_val_adj,
                   everything())
    }
    
    ## Select count percent preference---------------------
    if (p_o_c == "percent") {
        cell_type_out <- cell_type_out[, -grep("ount",
                                               colnames(cell_type_out))]
    }
    
    if (p_o_c == "count") {
        cell_type_out <- cell_type_out[, -grep("ercent",
                                               colnames(cell_type_out))]
    }
    
    
    if (drop_na_row == TRUE) {
        cell_type_out <- cell_type_out %>% drop_na()
    }
    
    ## Simplify stained list drop
    if (stained_gene_data == FALSE) {
        cell_type_out <- cell_type_out[, -grep("ected",
                                               colnames(cell_type_out))]
    }
    if (tested_protein_column == FALSE) {
        cell_type_out <- cell_type_out[, -grep("ested",
                                               colnames(cell_type_out))]
    }
    
    ## Final changes Make adjusted pvalues optional
    if (adjusted_pvals == FALSE) {
        cell_type_out <- cell_type_out[, -grep("dj", colnames(cell_type_out))]
    }
    
    
    
    
    
    ## Returning final table -------------
    return((cell_type_out))
}
