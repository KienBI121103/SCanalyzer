# Quality_Control.R

# Quality control for single - cell RNA - seq data

# Objective
# This script performs quality control on single-cell RNA-seq data using Seurat to remove low - quality cells

# Steps
# 1. Automation of quality control using ddqc
# 2. Remove of doublet cells using DoubletFinder
# 3. Quantify and remove cell free mRNAs from droplet based scRNA-seq data

# Loading required libraties
library(Seurat)
library(tidyverse)
library(ddqcR)
library(DoubletFinder)
library(SoupX)

# Loading objects requird for Quality Control
