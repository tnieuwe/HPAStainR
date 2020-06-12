# HPAStain.R

HPAStain.R is an R based package/Shiny app used to query the Human Protein Atlas for staining data. The purpose of this tool is to test if a list of proteins/genes is associated with a certain cell type in a HPA tested tissue. E.g. you have a list of protein coding genes from a differential expression single cell analysis and want to see if these proteins are associated with a known cell type. Instead of querying HPA multiple times you can load your list in HPAStainR which will return a ranked table of the cell types with the most protein staining.

This is my first time developing a package on Github, any and all feedback is appreciated!

## Getting Started

HPAStainR in its current form is a online Shiny app and R package. The package is not available on CRAN or Bioconductor, but that is our next step, until then feel free to download the binary from the github and install it.

Another way to use it is go to my shiny app website posted here and above:
https://32tim32.shinyapps.io/HPAStainR/ 

### Prerequisites

The R packages `tidyverse` and `shiny` are required for HPAStainR

The datasets required are below, but HPAStainR package has a function, `HPA_data_downlader()` that can download them for you either permanently or temporarily (`save_file` argument):

Normal Tissue: https://www.proteinatlas.org/download/normal_tissue.tsv.zip

Cancer Tissue: https://www.proteinatlas.org/download/pathology.tsv.zip

*Note*: these are large files but required to run HPAStainR

## Built With

* [Shiny](https://shiny.rstudio.com/) - For the shinyapp UI
* [Tidyverse](https://www.tidyverse.org/) - Data manipulation


## Authors

* **Tim Nieuwenhuis** - *Concept and Coder* - [32tim32](https://github.com/32tim32/)
* **Marc Halushka** - *Mentor* - [mhalushka](https://github.com/mhalushka)


## Acknowledgments

* Thank you Stephanie Yang for helpful comments - [syang93](https://github.com/syyang93/)
* HPA for data availability and data structure
* PurpleBooth for the README template - [PurpleBooth](https://gist.github.com/PurpleBooth/)
* Anyone else who has tested this tool
