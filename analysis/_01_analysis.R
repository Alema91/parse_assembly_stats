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
samples_ref <- read.table(paste0(path, "/data/samples_ref.txt"), header = F); colnames(samples_ref)<- c("id", "ref")
samples_id <- read.table(paste0(path, "/data/samples_id.txt"), header = F); colnames(samples_id)<- c("id")
fastq_path <- c("/srv/fastq_repo/MiSeq_GEN_317_20220713_MDFernanadez")


list_assembly<- list(0)
for (i in 1:length(samples_id)) {

    # Run, user, host and sequence
    name_run <- str_split(fastq_path, "/", simplify = T)[, 4]
    name_user <- str_split(posible_path, "_", simplify = T)[, 5]
    name_host <- tolower(str_split(posible_path, "_", simplify = T)[, 9])
    name_sequence <- samples_ref$ref[i]
    name_id <- samples_id$id[i]

    # columnas
    columnas <- "run\tuser\thost\tVirussequence\tsample\ttotalreads\treadshostR1\treadshost\t%readshost\tNon-host-reads\tContigs\tLargest_contig\t%Genome_fraction"
    name_columns <- as.vector(str_split(columnas, "\t", simplify = T))

    # totalreads
    json_fastp <- fromJSON(paste0("data/", name_sequence, "_20220726_viralrecon_mapping/fastp/", name_id, ".fastp.json"))
    value_totalreads <- json_fastp$summary[["after_filtering"]]$total_reads

    # readshostR1
    table_kraken <- read.table(paste0("data/", name_sequence, "_20220726_viralrecon_mapping/kraken2/", name_id, ".kraken2.report.txt"), sep = "\t")
    value_readhostr1 <- table_kraken$V2[table_kraken$V5 == 1]

    # readshosh
    value_readhost <- value_readhostr1 * 2

    # readshosh
    value_percreadhost <- table_kraken$V1[table_kraken$V5 == 1]

    # non host reads
    value_nonhostreads <- value_totalreads - value_readhost

    # % non host
    value_percnonhostreads <- round((value_readhost * 100) / value_totalreads, 2)

    # Contigs
    table_quast <- read.csv2("data/", name_sequence, "_20220726_viralrecon_mapping/assembly/spades/rnaviral/quast/transposed_report.tsv", skip = 0, sep = "\t", header = T)
    table_quast$id <- str_split(table_quast$Assembly, ".scaffolds", simplify = T)[, 1]
    table_ref_quast<- join(table_quast, samples_ref, by = "id")

    value_contigs <- table_quast$X..contigs[table_quast$id == name_id[i] | [table_quast$ref == name_ref[i]]
    value_lcontig <- table_quast$Largest.contig[table_quast$id == name_id[i]| [table_quast$ref == name_ref[i]]
    value_genomef <- table_quast$Genome.fraction....[table_quast$id == name_id[i]| [table_quast$ref == name_ref[i]]

    # Create table
    df_assembly <- data.frame(matrix(0, ncol = length(name_columns)))
    colnames(df_assembly) <- name_columns

    df_assembly[1, ] <- c(name_run, name_user, name_host, name_sequence, value_totalreads, value_readhostr1, value_readhost, value_percreadhost, value_nonhostreads, value_percnonhostreads, value_contigs, value_lcontig, value_genomef)

    list_assembly[[i]]<- df_assembly

}


