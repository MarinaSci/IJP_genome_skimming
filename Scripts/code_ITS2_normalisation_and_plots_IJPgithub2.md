---
title: "Code_ITS2_normalisation_plots_IJP"
author: "Marina Papaiakovou"
date: '2022-10-24'
output: html_document
---


## ITS2 analysis 
# Formatting of raw files and pulling ITS2 database
```{r, eval=FALSE}
# Genome skimming: code for normalisation of read data


- Downloaded complete ITS2 database from https://www.nemabiome.ca/its2-database.html on 5th Oct 2022
- Generated a file called "dada2.fasta" containing 11763 sequences
- there are many different species, but lots of duplication of some species.
- want to deduplicate this data so there is once representative sequence per species.
- difficult to know how to do this properly, ie, what is the "best/most representative" sequence? Decided to keep it simply and just choose the first sequence per species
```

```{r, eval=FALSE}
# counting unique sequence names
grep ">" dada2.fasta | sort | uniq | wc -l
#> 1671 unique species names


# create a new database, containing one fasta sequence per species
fastaq to_fasta -l0 dada2.fasta dada2_l0.fasta

grep ">" dada2_l0.fasta | sort | uniq | awk '{print $1}' |\
     while read -r NAME; do
          grep -A1 -m1 ${NAME} dada2_l0.fasta >> reference.fasta;
          done

# found that Marshallagia was causing false positives against bacterial reads, so removed it from the reference          
cat reference.fasta | sed '/Marshallagia/,+1d' >  reference.curated.fasta

```

```bash
#simple script to run each sample against ITS2 database using Vsearch.
#Requires changing the sample name and read files
```


# Vsearch for ITS2 mapping 



```{r, eval=FALSE}

sample_name=HA9
read1=HA9_trimmed_1.fq.gz
read2=HA9_trimmed_2.fq.gz

database=reference.curated.fasta

zcat $read1 > ${sample_name}_1.fastq
zcat $read2 > ${sample_name}_2.fastq

fastaq to_fasta ${sample_name}_1.fastq ${sample_name}_1.fasta
fastaq to_fasta ${sample_name}_2.fastq ${sample_name}_2.fasta

cat ${sample_name}_1.fasta r${sample_name}_2.fasta > ${sample_name}.fasta

rm ${sample_name}_1.fastq ${sample_name}_2.fastq ${sample_name}_1.fasta ${sample_name}_2.fasta

# run vsearch
bsub.py 10 --threads 10 vsearch_${sample_name} \
     "vsearch --usearch_global ${sample_name}.fasta \
     --db ${database} \
     --id 0.95 \
     --strand both \
     --blast6out ${sample_name}.vsearch_out \
     --minseqlength 100 \
     --mincols 100 \
     --maxaccepts 3 \
     --threads 7"

# sample_name=L1B
# database=reference.curated.fasta
#
#           bsub.py 10 --threads 10 vsearch_${sample_name} \
#                "vsearch --usearch_global ${sample_name}.fasta \
#                --db ${database} \
#                --id 0.95 \
#                --strand both \
#                --blast6out ${sample_name}.${database}.vsearch_out_min50 \
#                --minseqlength 50 \
#                --mincols 50 \
#                --maxaccepts 3 \
#                --threads 7"


```

## Human and livestock mitochondrial DNA analysis
```{r, eval=FALSE}
#R studio
# load libraries
library(tidyverse)

#Human mitogenome analysis 
# load data containing raw read counts, and fix the column names, and convert to a long, tidy format
hum_mito_raw_counts <- read.table("human_mito_rawcounts", sep="\t", header=F)

colnames(hum_mito_raw_counts) <- c("species", "start_pos", "genome_size_bp", "ET018", "ET103", "CAM1", "HA15", "HA16I", "HA9", "KA4E", "M104", "MA1G", "NDK113", "NDK63", "NDK92", "RA5C", "TKU102", "TKU23", "TKU25", "CAM2")

hum_mito_raw_counts_l <- pivot_longer(hum_mito_raw_counts, names_to = "sample_id", values_to = "raw_read_counts", cols=4:20)


# load data containing the total number of raw reads per sample
hum_samples_id_reads_n <- read.table("sample_id_raw_reads_n_human.txt", header=T, sep="\t")


# merged the dataframes
hum_mito_data <- full_join(hum_mito_raw_counts_l, hum_samples_id_reads_n, by="sample_id")


# normalise the data
hum_mito_data <- hum_mito_data %>% mutate(normalised = (raw_read_counts) / (raw_reads_n / 1000000) / (genome_size_bp / 1e6))


# remove counts lower than 10. Using this as a threshold to reduce the noise
hum_mito_data <- hum_mito_data %>% mutate(normalised = ifelse(normalised<10,0,normalised))


# make a plot
plot_hum_mito <- ggplot(hum_mito_data, aes(sample_id, species, fill=normalised)) +
     geom_tile(color = "gray") +
     theme_bw() +
     labs(x = "Species", y="Reads mapped / million reads / Mb", fill="Normalised\ncoverage") +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(face = "italic"))+
     scale_fill_viridis_c(option = "B", direction = -1, trans="log10",  na.value = "grey95", limits=c(7,50000)) +
     scale_y_discrete(limits=rev)

plot_hum_mito

## Livestock mitochondrial DNA analysis

liv_mito_raw_counts <- read.table("livestock_mito_rawcounts.txt", sep="\t", header=F)
colnames(liv_mito_raw_counts) <- c("species", "start_pos", "genome_size_bp", "L1B", "L2B", "L3B")

liv_mito_raw_counts_l <- pivot_longer(liv_mito_raw_counts, names_to = "sample_id", values_to = "raw_read_counts", cols=4:6)

liv_samples_id_reads_n <- read.table("sample_id_raw_reads_n_livestock.txt", header=T, sep="\t")

# merged the dataframes, first by "species" to add genome sizes, and then by "sample_id" to add the raw read counts
liv_mito_data <- full_join(liv_mito_raw_counts_l, liv_samples_id_reads_n, by="sample_id")

# normalise the data
liv_mito_data <- liv_mito_data %>% mutate(normalised = (raw_read_counts) / (raw_reads_n / 1000000) / (genome_size_bp / 1e6))
liv_mito_data <- liv_mito_data %>% mutate(normalised = ifelse(normalised<7,0,normalised))

plot_liv_mito <- ggplot(liv_mito_data, aes(sample_id, species, fill=normalised)) +
     geom_tile(color = "gray") +
     theme_bw() +
     labs(x = "Species", y="Reads mapped / million reads / Mb", fill="Normalised\ncoverage") +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(face = "italic"))+
     scale_fill_viridis_c(option = "B", direction = -1, trans="log10",  na.value = "grey95", limits=c(7,50000)) +
     scale_y_discrete(limits=rev)


plot_liv_mito


# join the human and livestock plots together. Due to the different sample sizes, using different widths to make the size of the samples approximately the same.
plot_hum_mito + plot_liv_mito + plot_layout(width=c(5.5,1), guides = "collect")

# save it
ggsave("figure2_mitochondrial_genome_normalised_cov.pdf", width=250, height=100, units="mm")

```


## Human whole genome analysis
```{r, eval=FALSE}
# read in the raw data per species, adding the species name as a new column
ascaris <- read.table("human_wholegenome/ascaris_genome_10kb.bed_all_stats", header=F)
ascaris$species <- "Ascaris lumbricoides"

anisakis <- read.table("human_wholegenome/anisakis_genome_10kb.bed_all_stats", header=F)
anisakis$species <- "Anisakis simplex"

ceylanicum <- read.table("human_wholegenome/ceylanicum_genome_10kb.bed_all_stats", header=F)
ceylanicum$species <- "Ancylostoma ceylanicum"

duodenale <- read.table("human_wholegenome/duodenale_genome_10kb.bed_all_stats", header=F)
duodenale$species <- "Ancylostoma duodenale"

echinococcus <- read.table("human_wholegenome/echinococcus_genome_10kb.bed_all_stats", header=F)
echinococcus$species <- "Echinococcus granulosus"

enterobius <- read.table("human_wholegenome/enterobius_genome_10kb.bed_all_stats", header=F)
enterobius$species <- "Enterobius vermicularis"

necator <- read.table("human_wholegenome/necator_genome_10kb.bed_all_stats", header=F)
necator$species <- "Necator americanus"

schistosoma <- read.table("human_wholegenome/schistosoma_genome_10kb.bed_all_stats", header=F)
schistosoma$species <- "Schistosoma mansoni"

strongyloides <- read.table("human_wholegenome/strongyloides_genome_10kb.bed_all_stats", header=F)
strongyloides$species <- "Strongyloides stercoralis"

taenia <- read.table("human_wholegenome/taenia_genome_10kb.bed_all_stats", header=F)
taenia$species <- "Taenia solium"

trichuris <- read.table("human_wholegenome/trichuris_genome_10kb.bed_all_stats", header=F)
trichuris$species <- "Trichuris trichiura"


# join the datasets together
human_genome_raw_data <- bind_rows(ascaris, anisakis, ceylanicum, duodenale, echinococcus, enterobius, necator, schistosoma, strongyloides, taenia, trichuris)

human_genome_raw_data <- human_genome_raw_data %>% group_by(species) %>% summarise(across(V4:V20, mean))

colnames(human_genome_raw_data) <- c("species", "ET018", "ET103", "CAM1", "HA15", "HA16I", "HA9", "KA4E", "M104", "MA1G", "NDK113", "NDK63", "NDK92", "RA5C", "TKU102", "TKU23", "TKU25", "CAM2")

human_genome_raw_data_l <- pivot_longer(human_genome_raw_data, names_to = "sample_id", values_to = "mean_read_counts", cols=2:18)

genome_sizes <- read.table("/Users/sd21/Desktop/marina/genome_size_wholegenome_human.txt", header=T, sep="\t")

samples_id_reads_n <- read.table("/Users/sd21/Desktop/marina/sample_id_raw_reads_n_human.txt", header=T, sep="\t")


# merged the dataframes, first by "species" to add genome sizes, and then by "sample_id" to add the raw read counts
data <- full_join(human_genome_raw_data_l, genome_sizes, by="species")
data <- full_join(data, samples_id_reads_n, by="sample_id")


# normalise the data
data <- data %>% mutate(normalised = (mean_read_counts * genome_size_bp/10000) / (raw_reads_n / 1000000) / genome_size_mb)


plot1 <- ggplot(data, aes(sample_id, normalised)) +
     geom_col(fill="#440154") +
     facet_wrap(species~., ncol = 3) +
     theme_bw() +
     labs(x = "Species", y="Reads mapped / million reads / Mb") +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.text = element_text(face = "italic"))

ggsave("figure3_human_genome_normalised_cov.pdf", width=170, height=170, units="mm")

```


## Livestock whole genome analysis

```{r, eval=FALSE}
#Rstudio
# load libraries
library(tidyverse)
library(patchwork)


# load each dataset, adding species name
canadensis <- read.table("livestock_wholegenome/canadensis_genome_10kb.bed_all_stats", header=F)
canadensis$species <- "Echinococcus canadensis"

gigantica <- read.table("livestock_wholegenome/gigantica_genome_10kb.bed_all_stats", header=F)
gigantica$species <- "Fasciola gigantica"

granulosus <- read.table("livestock_wholegenome/granulosus_genome_10kb.bed_all_stats", header=F)
granulosus$species <- "Echinococcus granulosus"

haemonchus <- read.table("livestock_wholegenome/haemonchus_genome_10kb.bed_all_stats", header=F)
haemonchus$species <- "Haemonchus contortus"

hepatica <- read.table("livestock_wholegenome/hepatica_genome_10kb.bed_all_stats", header=F)
hepatica$species <- "Fasciola hepatica"

oesophagostomum <- read.table("livestock_wholegenome/oesophagostomum_genome_10kb.bed_all_stats", header=F)
oesophagostomum$species <- "Oesophagostomum dentatum"

opisthorchis <- read.table("livestock_wholegenome/opisthorchis_genome_10kb.bed_all_stats", header=F)
opisthorchis$species <- "Opisthorchis viverrini"

saginata <- read.table("livestock_wholegenome/saginata_genome_10kb.bed_all_stats", header=F)
saginata$species <- "Taenia saginata"

solium <- read.table("livestock_wholegenome/solium_genome_10kb.bed_all_stats", header=F)
solium$species <- "Taenia solium"

suum <- read.table("livestock_wholegenome/suum_genome_10kb.bed_all_stats", header=F)
suum$species <- "Ascaris suum"

trichinella <- read.table("livestock_wholegenome/trichinella_genome_10kb.bed_all_stats", header=F)
trichinella$species <- "Trichinella spirallis"

trichuris <- read.table("livestock_wholegenome/trichuris_genome_10kb.bed_all_stats", header=F)
trichuris$species <- "Trichuris suis"


# join the datasets  
liv_genome_raw_data <- bind_rows(canadensis, gigantica, granulosus, haemonchus, hepatica, oesophagostomum, opisthorchis, saginata, solium, suum, trichinella, trichuris)


# Calculate the mean coverage per species
liv_genome_raw_data <- liv_genome_raw_data %>% group_by(species) %>% summarise(across(V4:V6, mean))


# fix names
colnames(liv_genome_raw_data) <- c("species", "L1B", "L2B", "L3B")


# convert to a tidy, long format
liv_genome_raw_data_l <- pivot_longer(liv_genome_raw_data, names_to = "sample_id", values_to = "mean_read_counts", cols=2:4)


# read in the genome size and raw read counts for each species
genome_sizes <- read.table("genome_size_wholegenome_livestock.txt", header=T, sep="\t")
samples_id_reads_n <- read.table("sample_id_raw_reads_n_livestock.txt", header=T, sep="\t")


# merged the dataframes, first by "species" to add genome sizes, and then by "sample_id" to add the raw read counts
data <- full_join(liv_genome_raw_data_l, genome_sizes, by="species")
data <- full_join(data, samples_id_reads_n, by="sample_id")


# normalise the data
data <- data %>% mutate(normalised = (mean_read_counts * genome_size_bp/10000) / (raw_reads_n / 1000000) / genome_size_mb)


# make a plot
plot2 <- ggplot(data, aes(sample_id, normalised)) +
     geom_col(fill="#440154") +
     facet_wrap(species~., ncol = 3) +
     theme_bw() +
     labs(x = "Species", y="Reads mapped / million reads / Mb") +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.text = element_text(face = "italic"))

plot2

ggsave("supfigure_livestock_genome_normalised_cov.pdf", width=170, height=170, units="mm")
```


# Testing Enterobius contamination
- we found that small contigs contained higher coverage, and when blasted, presented as bacterial hits
- was to see if the true signal improves when removing contigs smaller than 10kb
- aim is to compare data before and after filtering

```{r, eval=FALSE}
#R studio 
library(tidyverse)

# read in the enterobius data
enterobius <- read.table("human_wholegenome/enterobius_genome_10kb.bed_all_stats", header=F)
enterobius$species <- "Enterobius vermicularis"


# keep all rows above row containing "4528_enterobius_vermicularis", which is th last scaffold above 10 kb
enterobius_filtered <- enterobius %>% filter(cumsum(grepl("4528_enterobius_vermicularis",V1,fixed=TRUE))<1)
enterobius_filtered$species <- "Enterobius vermicularis - above 10 kb"


# join the two datasets
enterobius_genome_raw_data <- bind_rows(enterobius, enterobius_filtered)


# calculate mean per group
enterobius_genome_raw_data <- enterobius_genome_raw_data %>% group_by(species) %>% summarise(across(V4:V20, mean))


# fix column names as sample IDs
colnames(enterobius_genome_raw_data) <- c("species", "ET018", "ET103", "CAM1", "HA15", "HA16I", "HA9", "KA4E", "M104", "MA1G", "NDK113", "NDK63", "NDK92", "RA5C", "TKU102", "TKU23", "TKU25", "CAM2")


# convert to long, tidy format
enterobius_genome_raw_data_l <- pivot_longer(enterobius_genome_raw_data, names_to = "sample_id", values_to = "mean_read_counts", cols=2:18)


# import raw read counts per sample
samples_id_reads_n <- read.table("/Users/sd21/Desktop/marina/sample_id_raw_reads_n_human.txt", header=T, sep="\t")


# merged the dataframes, by "sample_id" to add the raw read counts
data <- full_join(enterobius_genome_raw_data_l, samples_id_reads_n, by="sample_id")


# add the number of 10 kb genome windows in each group
data <- data %>% mutate(genome_windows = ifelse(species=="Enterobius vermicularis", 28858, ifelse(species=="Enterobius vermicularis - above 10 kb", 13552, 0)))


# normalise the data
data <- data %>% mutate(normalised = (mean_read_counts * genome_windows) / (raw_reads_n / 1000000) / (genome_windows*10000 / 1e6))


# make a plot
plot_enterobious <- ggplot(data, aes(sample_id, normalised)) +
     geom_col(fill="#440154") +
     facet_wrap(species~., ncol = 3, scales = "free") +
     theme_bw() +
     labs(x = "Species", y="Reads mapped / million reads / Mb") +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.text = element_text(face = "italic"))

plot_enterobious

# save it
ggsave("SupplementaryFigure_enteriobius.pdf", width=170, height=100, units="mm")

```

# Diagnostics comparison plot 
Comparison between genome skimming and PCR plot, using 'pie charts' per sample, per helminth  tested
```{r}
ggplot(data, aes(x="", y="", group=dx, colour=dx, fill=dx, alpha=dx_result)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  facet_grid(helminth~sample_ID) +
  guides(color = "none") +
  labs(fill="Diagnostic method", y="Sample ID", x="Species", alpha="Diagnostic result", ylab(""), xlab("")) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(), strip.text.y = element_text(angle = 360, face="italic"), strip.text.x = element_text(angle = 90))+ #this will flip the names of the y, x axis 
  
  theme(legend.position="bottom", legend.box = "horizontal")+
  scale_color_viridis_d(option = "cividis")+
  scale_fill_viridis_d(option = "cividis")

ggsave("Figure_5_dx_plots.pdf", useDingbats = FALSE, width = 210, height = 110, units = "mm")

```

