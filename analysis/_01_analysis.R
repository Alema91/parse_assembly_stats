#!/usr/bin/env Rscript

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(plyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(readxl, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(forcats, quietly = TRUE, warn.conflicts = FALSE)
library(reshape2, quietly = TRUE, warn.conflicts = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(magick, quietly = TRUE, warn.conflicts = FALSE)
library(kableExtra, quietly = TRUE, warn.conflicts = FALSE)
library(knitr, quietly = TRUE, warn.conflicts = FALSE)
library(magrittr, quietly = TRUE, warn.conflicts = FALSE)
library(ggpubr, quietly = TRUE, warn.conflicts = FALSE)

library(jsonlite, quietly = TRUE, warn.conflicts = FALSE)

### Not in o in ----

`%notin%` <- Negate(`%in%`)

################################################
################################################
## DATA          ###############################
################################################
################################################

# PATHS
path <- getwd()
posible_path <- c("/data/bi/scratch_tmp/bi/SRVCNM716_20220722_GENOMEEV15_mdfernandez_S/ANALYSIS/20220726_ANALYSIS01_METAGENOMIC_HUMAN")
samples_ref <- read.table(paste0(path, "/data/samples_ref.txt"), header = F)
samples_id <- read.table(paste0(path, "/data/samples_id.txt"), header = F)

# Run, user and host
name_run <- str_split(posible_path, "/", simplify = T)[, 6]
name_user <- str_split(posible_path, "_", simplify = T)[, 5]
name_host <- tolower(str_split(posible_path, "_", simplify = T)[, 9])

# columns names
columnas <- "run\tuser\thost\tVirussequence\tsample\ttotalreads\treadshostR1\treadshost\t%readshost\tNon-host-reads\t%Non-host-reads\tContigs\tLargest_contig\t%Genome_fraction"
name_columns <- as.vector(str_split(columnas, "\t", simplify = T))

# Total reads
json_fastp <- fromJSON("data/EVD68_20220726_viralrecon_mapping/fastp/PTA102.fastp.json")
value_totalreads <- json_fastp$summary[["after_filtering"]]$total_reads

# readshostR1
table_kraken <- read.table("data/EVD68_20220726_viralrecon_mapping/kraken2/PTA102.kraken2.report.txt", sep = "\t")
value_readhostr1 <- table_kraken$V2[table_kraken$V5 == 1]

# readshosh
value_readhost <- value_readhostr1 * 2

# readshosh
value_percreadhost <- table_kraken$V1[table_kraken$V5 == 1]
