library(tidyverse)
library(gridExtra)
library(reshape2)
library(ggpubr)
library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(shinyjs)

aa_list <- c("A","C", "D", "E", "F", "G", "H", "I", "K" ,"L", "M", "N" ,"P", "Q" ,"R" ,"S" ,"T", "V" ,"W" ,"Y")

main_df <- read_csv("../data/main_df.csv", na=character())
