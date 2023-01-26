# MCLClusterJCI
Code and input file to regenerate the cluster result of publication: 
https://www.jci.org/articles/view/153283/ga
The original cluster code is from https://github.com/broadinstitute/DLBCL_Nat_Med_April_2018

Shiny App of project MCL samples into 4 clusters identified in publication: 


Genetic lesion reference matrix also included samples from https://ashpublications.org/blood/article/136/12/1419/461130/Genomic-and-epigenomic-insights-into-the-origin


Method1: Load this link into your browser.
Shiny App link: https://meilingjin1002.shinyapps.io/MCLClusterJCI/

Method2: Run from github
Open your Rstudio

library(shiny)

runGitHub(repo='JCI2022', username='COHWangLab', ref = 'main')



The example input is data/example_input.txt.
Please prepare the input file identical to the data/example_input.txt with rows are genetic lesions, and columns are sample ids. 
1. Upload your input file.
2. Select parameters.

