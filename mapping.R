#IJP analysis of genome skimming data, using illumina short reads mapped directly to master reference files ----

#Databases ----
#Human reference files
#mitogenome database includes: 
AP017684_Enterobius_vermicularis_mitochondrial_DNA_complete
NC_002545_Schistosoma_mansoni_mitochondrion_complete_genome
NC_003415_Ancylostoma_duodenale_mitochondrion_complete_genome
NC_003416_Necator_americanus_mitochondrion_complete_genome
NC_004022_Taenia_solium_mitochondrion_complete_genome
NC_007934_Anisakia_simplex_mitochondrion_complete_genome
NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome
NC_017750_Trichuris_trichiura_mitochondrion_complete_genome
NC_028624_Strongyloides_stercoralis_isolate_PV001_mitochondrion
NC_035142_Ancylostoma_ceylanicum_mitochondrion_complete_genome
NC_044548_Echinococcus_granulosus_mitochondrion_complete_genome

#Human reference files -----
#whole genomes
#all from wormbase
Ascaris_lumbricoides
Necator_americanus
Trichuris_trichiura
Strongyloides_stercoralis
Ancylostoma_ceylanicum
Taenia_solium
Schistosoma_mansoni
Echinococcus_granulosus
Enterobius_vermicularis
Anisakis_simplex

#Livestock reference files ----
#livestock mitogenome files (NCBI):
Schistosoma_curassoni
Schistosoma_japonicum 
Ascaris_suum
Trichuris_suis
Fasciola_hepatica
Trichinella_spirallis
Taenia_solium
Taenia_asiatica
Haemonchus_contortus
Trichostrongylus_vitrinus
Oesophagostomum_dentatum
Trichostrongylus_axei
Fasciola_gigantica

#livestock whole genome reference files, all from wormbase

Echinococcus_granulosus
Trichinella_spirallis
Haemonchus_contortus
Fasciola_hepatica
Fasciola_gigantica
Ascaris_suum
Trichuris_suis
Oesophagostomum_dentatum
Opisthorchis_viverrini
Taenia_saginata
Echinococcus_canadensis

#livestock ribosomal (because no full genome is available)
Trichostrongylus_spp
Dicrocoelium_spp
Camelostrongylus 
Glanothostoma_spp 
Hyostrongylus_spp 
Mecistocircus spp 
Cysticercus_spp
Bertiella_spp
Oesophagostum_spp
Orloffia_spp 

#Renaming the full genome contig files to #number_of_contig_species_name
#Did not figure out how to do this in a loop so I did them individually 
#make sure you are inside a conda -env that allows you to use fastaq 
fastaq enumerate_names refence_file.fa reference_file_renamed.fa --suffix _speciesname
#for a loop
for i in *.fa; do fastaq enumerate_names $i ${i}sta --suffix ${i}; done #i don't think this works ok 


#once you have renamed all the > contigs for all the files, put them into a master file 
cat *_renamed.fa >> human_parasite_genomes.fasta 

#going back to my original script for variant calling: 
#SCRIPT for mass mapping ----
-----------------------------------------------------------------
#!/bin/env
set -e #to hide any errors
#where the fasta file is and index that
reference=your_ref.fasta
#set prefix
prefix=project_name #or experiment name

bwa index $reference

for fq in *_1.fq.gz
do
      echo "working with $fq"
      base=$(basename -a $fq | sed 's/_1.fq.gz//g')
      echo ""
      
      echo "genome file is: $reference"
      echo "base name is: $base"
      
      base_1=${base}_1.fq.gz
      base_2=${base}_2.fq.gz
      
      echo "base_1 file is: $base_1"
      echo "base_2 file is: $base_2"
      
      echo ""
      echo "I will now run the following command:"
      echo "bwa mem ${reference} ${base_1} ${base_2} > ${prefix}_${base}.sam"

      bwa mem ${reference} <(zcat ${base_1}) <(zcat ${base_2}) > ${prefix}_${base}.sam  #seems to have worked
      #tried some stringent mapping, bwa mem -B 40 -O 60 -E 10 -L100 human_mito_ref.fasta <(zcat TKU102_trimmed_1.fq.gz) <(zcat TKU102_trimmed_2.fq.gz) > TKU102_stringent.sam
      #but decreased those that properly paired by n=700 or so 
      
      #STILL NEEDS SOME WORK BELOW - IT IS NOT FINAL 
      #samtools view -q 30 -F 4 -S -b ${prefix}_${base}.sam > ${prefix}_${base}.bam #think about  adding -q 30 - q-10 changed the mapping by 
      #rm ${prefix}_${base}.sam
      #samtools sort -o ${prefix}_${base}_sorted.bam ${prefix}_${base}.bam
      #rm ${prefix}_${base}.bam
      #maybe samtools view here flagstat
      #samtools index ${prefix}_${base}_sorted.bam #consider add prefix here
      #samtools idxstats ${prefix}_${base}_sorted.bam >> ${prefix}_${base}.txt
      #or samtools coverage etc - only in most up to date versions of samtools 
      #samtools coverage -l 135 -q 40 -Q 40 -d 400 -b TKU102_list -m
done

----------------------------------------------------------------------
#assumes you have indexed your reference and have run the above cod, before the comments 
  
  
#01.ALIGN ----
bwa mem ${reference} <(zcat ${base_1}) <(zcat ${base_2}) > ${prefix}_${base}.sam  

  
#SAM FILES and how to handle them 
#what is the output of a sam file 
#if you are getting a lot of 0 * 0 * then maybe you are looking at unmapped reads!!! 
ST-E00523:583:HYCNYCCXY:4:1101:8866:1977        83      NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome(CHROM)    13181   60      134M (CIGAR)    =       12902   -413    
TTTTTTTTTGTTCGTTCTAATGTAATTTTAATATATATGTGTTTTGAGTTGTCCATGTTTCCTATTTTGGTTATAATTCTTGGCTATCGTTCTTAGATTGAGAAAATTAATTCTTCTTATTATTTAATTTTTTA  
FJJFFJAF--A---F7-J7-FFA7-FF<A--A-A<<7FAFJFJ<F<FAFJAA-A-FJJJ<-AF-JJAJFJJFFAJJFFAAF7<-FJFFJJF---FJJJJJJJFAJJFFJJF-JJJJJJJFJJFJJAJJJJJJJF  NM:i:11 MD:Z:13C12C2G2G7A12T5C23T3G4C0C40       MC:Z:73M60S     AS:i:79 XS:i:0  
http://bio-bwa.sourceforge.net/bwa.shtml#4
#	MD: Mismatching positions/bases
#	AS: Alignment score
#XS:Suboptimal alignment score
1.Read Name
2.SAM flag
3.chromosome (if read is has no alignment, there will be a “*” here)
4.position (1-based index, “left end of read”)
5.MAPQ (mapping quality — describes the uniqueness of the alignment, 0=non-unique, >10 probably unique)
6.CIGAR string (describes the position of insertions/deletions/matches in the alignment, encodes splice junctions, for example)
7.Name of mate (mate pair information for paired-end sequencing, often “=”)
8.Position of mate (mate pair information)
9.Template length (always zero for me)
10.Read Sequence
11.Read Quality
12.Program specific Flags (i.e. AS is an alignment score, NH is a number of reported alignments that contains the query in the current record)

#the CIGAR is something that apparently is super helpful for my data! Because I am getting partial mapping (e.g., 30 bp out of 130 bp read that maps something non specific)
#using something like samclip below, you get rid of some very hard alignments. For example, 92M1D43M  means that 92 mapped on one side, then there was a deletion and then 43 mapped on the other side 
#this was retained by the program but hard clipping was not. 


#02.Remove unmapped reads ----
#get rid of unmapped reads
samtools view -q 30 -F 4 -S -h ${prefix}_${base}.sam > ${prefix}_${base}_onlymapped.sam
#testing 

#03.Filter reads by length ----
samtools view -h ${prefix}_${base}_onlymapped.sam | awk 'length($10) > 130 || $1 ~ /^@/' > ${prefix}_${base}_onlymapped_filtered.sam

#samclip for CIGAR ----
#Used conda to install samclip, that gets rid of super hard clipped alignments in a sam file... 
conda install -c bioconda -c conda-forge samclip

#04.Index FASTA for samclip ----
#samclip needs a ref file that has been indexed differently:
samtools faidx human_mito_ref.fasta


#output of /fai file: 
#NAME: name of ref 
#LENGTH: length of seq in bases
#OFFSET: offset within the FASTA file of this sequence first base
#LINEBASES: number of bases in each line
#LINEWIDTH: in bytes

  
#run samclip 
samclip --ref ${reference}.fai ${prefix}_${base}_onlymapped_filtered.sam > ${prefix}_${base}_samclip.sam

#05.SAM TO BAM ----
#convert
samtools view -S -b ${prefix}_${base}_samclip.sam > ${prefix}_${base}.bam
#sort
samtools sort -o ${prefix}_${base}_sorted.bam ${prefix}_${base}.bam

#remove duplicates
sambamba view -t 12 -h -f bam -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" ${prefix}_${base}_sorted.bam -o ${prefix}_${base}_uniq.bam

#index the bam 
samtools index ${prefix}_${base}_sorted.bam

#06.Calculate depth for coverage plots ----

#calculate depth per base
samtools depth ${prefix}_${base}_sorted.bam > deduped_${prefix}_${base}.coverage

#07.Visualise coverage (alternatively)
samtools coverage -b .bam -m
#bam required as a separate file/list, however it does not do per sample coverage, that's why I did the above 
#-m will give you stdout 
# > .tsv will give you a table delimited file that you can merge/plot later 
#(steves_awesome_environment) marip3@franklin:~/mbl_genome_skimming/ijp_human/test$ less TKU102_samclip_onlymapped_filtered_sorted.tsv 

#output of samtools coverage  
#rname  startpos        endpos  numreads        covbases        coverage        meandepth       meanbaseq       meanmapq
AP017684_Enterobius_vermicularis_mitochondrial_DNA_complete     1       14003   2       39      0.278512        0.00278512      34.9    40
NC_003415_Ancylostoma_duodenale_mitochondrion_complete_genome   1       13721   6       23      0.167626        0.00845419      38.7    42
NC_004022_Taenia_solium_mitochondrion_complete_genome   1       13709   5       41      0.299074        0.0070027       37.9    41.2
NC_007934_Anisakia_simplex_mitochondrion_complete_genome        1       13916   8       39      0.280253        0.0109945       30.4    40.8
NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome    1       14281   19452   14222   99.5869 182.558 37.5    60
NC_035142_Ancylostoma_ceylanicum_mitochondrion_complete_genome  1       13660   3       290     2.12299 0.0212299       39.2    53.3
NC_002545Schistosoma_mansoni_mitochondrion_complete_genome      1       14415   0       0       0       0       0       0
NC_003416_Necator_americanus_mitochondrion_complete_genome      1       13605   0       0       0       0       0       0
NC_017750_Trichuris_trichiura_mitochondrion_complete_genome     1       14046   0       0       0       0       0       0
NC_028624_Strongyloides_stercoralis_isolate_PV001_mitochondrion 1       13751   0       0       0       0       0       0
NC_044548_Echinococcus_granulosus_mitochondrion_complete_genome 1       17675   0       0       0       0       0       0
TKU102_samclip_onlymapped_filtered_sorted.tsv (END)



#I still get some in ceylanicum so I need to see what these are 
#update - only get a couple on ceylanicum that probably would be a matter of ID match, maybe if I select them as unique 

#Need to figure out if there are still things mapping to ceylanicum for example even partially that way - 
#update: no, 1-2 reads are left, so the pipeline is good enough for now 



#USEFUL TIPS TO PRODUCE COVERAGE PLOTS/STATS ------

#AWK for results with > 10 reads 
cat yourfile_trimmed_uniq.txt | awk -F "\t" '{ if($3 >10) { print } }' > yourfile_trimmed_uniq_onlyresultswithreads.txt

#FASTA.FAI to BED FILE -----
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' yourfile.fasta.fai > yourfile.bed

#BEDTOOLS MULTICOV for bam stats ----
bedtools multicov -bams bam1 bam2 bam3 -bed reference.bed > yourfile_stats.txt
bedtools multicov -bams *_uniq.bam -bed human_mito.bed > stats_mito_all.txt
#if you run the wild card, then the bam files should have been taken alphabetically by the software 
#to verify 
ls -1 *_uniq.bam

#From Steve: 
#Each row in the output of stats.mito)all.txt is a contig - in this case a whole mtDNA genome from each speciess.
#First column is the name, second cant remember, third is sequence length, and fourth onward is the bases mapped for each BAM file
#I don't get a header, but I can use this in a heamap 

#Then you can use something powerful like "datamash" to generate the summary stats - super useful tool. 
https://www.gnu.org/software/datamash/

#Remove spaces from headers in references ----
conda install -c agbiome bbtools
reformat.sh in=livestock_ribo_ref.fasta out=livestock_ribo_ref_fixed.fasta addunderscore

#CUT OUT .bed file per species with contig name and column 1 and 3
cat yourfile.fasta.fai | head | grep -Fvw "yourparasite" | cut -f 1,2 > species.genome

#ISOLATE CONTIGS FOR ALL SPECIES IN A MULTIREFERENCE FILE, AND KEEP COLUMN 1 AND 2 ----
cut -f 1 livestock_parasite_full_genomes.fasta.fai | cut -f 2 -d "_" | sort | uniq | while read CHR; do grep $CHR livestock_parasite_full_genomes.fasta.fai | cut -f 1-2 > ${CHR}_genome; done
#it does not like it for the human genome stuff

#ceylanicum and duodenale ended up togerher... 
cut -f 1 ancylostoma_genome | cut -d'_' -f3 | sort | uniq | while read CHR; do grep $CHR ancylostoma_genome | cut -f 1-2 > ${CHR}_genome; done
#WORKED !!!!

#remove the .fa from the end of every line in column 1 
cat filename.fasta | sed 's/.fa//' > copy_edited

#BEDTOOLS makewindows ----
#from individual bed files, species-specific genomes from fasta.fai so I can do the multi cov 

bedtools makewindows -g duodenale_genome -w 10000 > 

#BEDTOOLS makewindows ----
#loop
for i in *_genome; do bedtools makewindows -g $i -w 10000 > ${i}_10kb.bed; done

#BEDTOOLS multicov ----
#per individual species genome, do multicov for all samples.... 
bedtools multicov -bed duodenale_genome_10kb.bed -bams *_uniq.bam > duodenale_all_stats

#loop bedtools multicov ----
for i in *_10kb.bed; do bedtools multicov -bed $i -bams *_uniq.bam > ${i}_all_stats; done
#finished for livestock too

#Normalise reads ----
#I would do the following:

#(reads mapped / total reads) *1,000,000
#where total reads is reads per sample
















#notes ----
#make sure you give it a list, even if you only have one file in there 
samtools coverage -l 100 -d 1000000 -b bamlist -m

#if you don't want the output in the screen, then do 
samtools coverage -l 100 -d 1000000 -b bamlist -o nameoffile_mitocoverage.tsv

#experimented with filters but I needed to sort out the filters in the sam files first (see above)
samtools coverage -l 135 -q 40 -Q 40 -d 400 -b TKU102_list -m

#random stuff I tried with samtools and bam files ----
samtools view -q 10 -F 1284 -f 0x02 -b TKU102_stringent_sorted.bam > TKU102_properlypaired.bam

samtools view -h TKU102_properlypaired.bam | awk 'length($10) > 130 || $1 ~ /^@/' | samtools view -bS - > TKU102_properlypaired_filtered.bam

samtools view TKU102_properlypaired_filtered.bam | perl -lane 'print if $F[5] =~ /^90M$/;' > TKU102_filtered100id.bam

java -jar dist/bioalcidaejdk.jar -e 'stream().map(R->R.getReadUnmappedFlag()?0:(int)(100.0*(R.getReadLength()-R.getIntegerAttribute("NM"))/(double)R.getReadLength())).collect(Collectors.groupingBy(Function.identity(), Collectors.counting())).forEach((K,V)->{println(K+"\t"+V);});' TKU102_properlypaired_filtered.bam 

