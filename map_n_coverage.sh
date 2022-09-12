#!/usr/env/bin bash
set -e #to hide any errors

#where my script is

#where the fasta file is and index that
reference=human_mito_ref.fasta
#set prefix
prefix=ijp_human_mito #or experiment name

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
        #i checked this and it works on its own
        bwa mem ${reference} <(zcat ${base_1}) <(zcat ${base_2}) > ${prefix}_${base}.sam  #seems to have worked

        #get rid of unmapped reads
        samtools view -q 30 -F 4 -S -h ${prefix}_${base}.sam > ${prefix}_${base}_onlymapped.sam #keep the header -h

        #filter them by length
        samtools view -h ${prefix}_${base}_onlymapped.sam | awk 'length($10) > 130 || $1 ~ /^@/' > ${prefix}_${base}_onlymapped_filtered.sam

        #if not done so already, index the fasta to be used with samclip

        samtools faidx $reference
        #then  filter them by hard clipping
        samclip --ref ${reference}.fai ${prefix}_${base}_onlymapped_filtered.sam > ${prefix}_${base}_samclip.sam
        #convert to bam
        samtools view -S -b ${prefix}_${base}_samclip.sam > ${prefix}_${base}.bam
        #sort
        samtools sort -o ${prefix}_${base}_sorted.bam ${prefix}_${base}.bam
        #remove duplicates

        sambamba view -t 12 -h -f bam -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" ${prefix}_${base}_sorted.bam -o ${prefix}_${base}_uniq.bam

        samtools index ${prefix}_${base}_uniq.bam

        samtools idxstats ${prefix}_${base}_uniq.bam > ${prefix}_${base}_uniq.txt

        #calculate depth per base
        samtools depth ${prefix}_${base}_uniq.bam > deduped_${prefix}_${base}.coverage

        #I know that ascaris worked best, however, I will need to adapt this for multiple chromosomes????
        #awk '$1 == "NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome" {print $0}' deduped_TKU102.coverage > ascaris_${prefix}_${base}.coverage
        #take table .coverage and plot it in R

        done


#ignore
#ls -1 *_trimmed_sorted.bam > sample_list

#samtools coverage -b sample_list > sample.tsv
