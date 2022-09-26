# parse_assembly_stats
Parsear los valores del output de viralrecon para obtener la tabla de assembly stats. Esta tabla consistiría en las siguientes columnas:

1. Run
2. User
3. Host
4. Sequence
5. Sample
6. Total reads
7. Read host R1
8. Read host
9. % Read host
10. Non host reads
11. % Non host reads
12. Contigs
13. Largest contig
14. % Genome fraction

# Modo de uso

```
# conda activate R-4.1.2

cp create_assembly_stats.R *_ANALYSIS01_METAGENOMIC_HUMAN
Rscript create_assembly_stats.R

```