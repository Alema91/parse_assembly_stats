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
library(writexl, quietly = TRUE, warn.conflicts = FALSE)

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
colnames(samples_ref) <- c("id", "ref")
fastq_path <- c("/srv/fastq_repo/MiSeq_GEN_317_20220713_MDFernanadez")

# columnas
columnas <- "run\tuser\thost\tVirussequence\tsample\ttotalreads\treadshostR1\treadshost\t%readshost\tNon-host-reads\t%Non-host-reads\tContigs\tLargest_contig\t%Genome_fraction"
name_columns <- as.vector(str_split(columnas, "\t", simplify = T))

list_assembly <- list(0)
for (i in 1:nrow(samples_ref)) {

    # Run, user, host and sequence
    name_run <- str_split(fastq_path, "/", simplify = T)[, 4]
    name_user <- str_split(posible_path, "_", simplify = T)[, 5]
    name_host <- tolower(str_split(posible_path, "_", simplify = T)[, 9])
    name_sequence <- samples_ref$ref[i]
    name_id <- samples_ref$id[i]

    # totalreads
    json_fastp <- fromJSON(paste0("data/", name_sequence, "_20220726_viralrecon_mapping/fastp/", name_id, ".fastp.json"))
    value_totalreads <- json_fastp$summary[["after_filtering"]]$total_reads

    # readshostR1
    table_kraken <- read.table(paste0("data/", name_sequence, "_20220726_viralrecon_mapping/kraken2/", name_id, ".kraken2.report.txt"), sep = "\t")
    value_readhostr1 <- table_kraken$V2[table_kraken$V5 == 1]

    # readshosh
    value_readhost <- value_readhostr1 * 2

    # readshost
    value_percreadhost <- table_kraken$V1[table_kraken$V5 == 1]

    # non host reads
    value_nonhostreads <- value_totalreads - value_readhost

    # % non host
    value_percnonhostreads <- round((value_readhost * 100) / value_totalreads, 2)

    # Contigs
    table_quast <- read.csv2(paste0("data/", name_sequence, "_20220726_viralrecon_mapping/assembly/spades/rnaviral/quast/transposed_report.tsv"), skip = 0, sep = "\t", header = T)
    table_quast$id <- str_split(table_quast$Assembly, ".scaffolds", simplify = T)[, 1]
    table_ref_quast <- join(table_quast, samples_ref, by = "id")

    value_contigs <- as.numeric(table_ref_quast$X..contigs[table_ref_quast$id == name_id & table_ref_quast$ref == name_sequence])
    value_lcontig <- as.numeric(table_ref_quast$Largest.contig[table_ref_quast$id == name_id & table_ref_quast$ref == name_sequence])
    value_genomef <- as.numeric(table_ref_quast$Genome.fraction....[table_ref_quast$id == name_id & table_ref_quast$ref == name_sequence])

    # empty values
    if (length(value_contigs) == 0) {
        value_contigs <- NA
        value_lcontig <- NA
        value_genomef <- NA
    }

    # Create table
    list_assembly[[i]] <- c(name_run, name_user, name_host, name_sequence, name_id, value_totalreads, value_readhostr1, value_readhost, value_percreadhost, value_nonhostreads, value_percnonhostreads, value_contigs, value_lcontig, value_genomef)
}

df_final <- as.data.frame(do.call("rbind", list_assembly))
colnames(df_final) <- name_columns

# Write table
write.table(df_final, "results/assembly_stats.csv", row.names = F, col.names = T, sep = "\t", quote = F)
