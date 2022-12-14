# MCLClusterJCI
Shiny App of project MCL samples into 4 clusters identified in publication: 
https://www.jci.org/articles/view/153283/ga



Method1: Load this link into your browser.
Shiny App link: https://meilingjin1002.shinyapps.io/MCLClusterJCI/

Method2: Run from github
Open your Rstudio
library(shiny)
runGitHub(repo='JCI2022', username='COHWangLab', ref = 'main')


The example input is data/example_input.txt.
Please prepare the input file identical to the data/example_input.txt with rows are genetic lesions, columns are sample ids. 
1. Upload your input file.
2. Select parameters.

