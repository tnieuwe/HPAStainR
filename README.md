# HPAStain.R

HPAStain.R is an R based tool used to query the Human Protein Atlas for staining data. The purpose of this tool is to test if a list of proteins is associated with a certain cell type in a tissue. E.g. you have a list of protein coding genes from a differential expression single cell analysis and want to see if these proteins are associated with a known cell type. Instead of querying HPA multiple times you can load your list in HPAStain.R which will return a ranked table of the cell types with the most protein staining.

This is my first time developing a tool on github, any and all feedback is appreciated!

## Getting Started

HPAStain.R in its current form is a shiny app an an R function. As it is a single function have not turned it into a package and so to use the function you can just copy or clone the function from the R code (working on streamlining this)

Another way to use it is go to my shiny app website posted here and above:
https://32tim32.shinyapps.io/HPAStainR/ 

### Prerequisites

To run the shinyapp locally you will need the shiny library otherwise you only need tidyverse and the two datasets below:

Normal Tissue: https://www.proteinatlas.org/download/normal_tissue.tsv.zip
Cancer Tissue: https://www.proteinatlas.org/download/pathology.tsv.zip

Note these are large files but required to run HPAStain.R, you could also make a pipeline with the package hpar if you are inclined to working completely in R. 


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
