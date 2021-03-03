#!/bin/bash
#PBS -l select=1:ncpus=10:mem=100gb:scratch_local=100gb
#PBS -l walltime=20:20:00
#PBS -N Tags_barcode01

###part0 - run scripts below, only ONCE
#prepare_ref.sh
#index_fast5.sh
#prepare_control_genome.sh


module add python-2.7.10-gcc
module add conda-modules-py37
conda activate /storage/brno2/home/nikolato/.conda/envs/test

USER=nik
DATE=$(date +'%m_%d_%Y')
ORGANISM=SS14
SAMPLE=BC01 ###SAMPLE DIR NAME; including fastq files
TEST_UNIT=ORF ### "tags" or "ORF"
CONTROL_GENOME=/storage/brno3-cerit/home/nikolato/ngs/data/TPE_tags/control_regions/$ORGANISM/$TEST_UNIT/OMP_ORF_reference_sequences_SS14  
FASTQ=/storage/brno3-cerit/home/nikolato/ngs/data/TPE_tags/Pilot/fastq/demultiplexed/$SAMPLE
RESULTS=$FASTQ/results/$USER/$TEST_UNIT/$DATE/$SAMPLE.pipeline4
REFERENCE=/storage/brno3-cerit/home/nikolato/ngs/references/references/treponema/ss14/CP004011.fasta
REF_DIR=/storage/brno3-cerit/home/nikolato/ngs/data/TPE_tags/reference/$ORGANISM
BED_FULL_GENOME=/storage/brno3-cerit/home/nikolato/ngs/data/TPE_tags/tags/$ORGANISM/full_genome.bed
BED=/storage/brno3-cerit/home/nikolato/ngs/data/TPE_tags/tags/$ORGANISM/Tagy_Nik_PCR_end_2020_08_26.bed
TAG_DIR=/storage/brno3-cerit/home/nikolato/ngs/data/TPE_tags/tags/$ORGANISM
INSERT_SIZE_FILE=/storage/brno3-cerit/home/nikolato/ngs/data/TPE_tags/tags/insert_size/$ORGANISM/Tagy_Nik_PCR_end_2020_08_26.insert_sizeo3-cerit/home/nikolato/ngs/data/TPE_tags/control_regions/$ORGANISM/ORF/OMP_ORF_reference_sequences_SS14

ALIGNMENT_MASK=$RESULTS/alignment_masked_ref
TAG_DIR_MASKED=$TAG_DIR/tags_masked
LONG_READS=$RESULTS/extracted_long_reads
ASSEMBLY=$RESULTS/assembly_masked_ref
MAPPING=$RESULTS/mapping_masked_ref
RAW_REF_CONTIGS=$RESULTS/raw_ref_contings_masked_ref
VARIANTS=$RESULTS/variants_masked_ref
POLISHED_GENOME=$RESULTS/polished_genome_masked_ref
ALL_READS=$FASTQ/merged.fastq
RACON=$RESULTS/racon/$SAMPLE
MEDAKA=$RESULTS/medaka/$SAMPLE
TEST_RACON=$RACON/performance_test
TEST_MEDAKA=$MEDAKA/performance_test

GENOME_SIZE="6k"


#########MM filtering#########
SAMTOOLS=samtools
BAMUTILS=/storage/brno3-cerit/home/nikolato/ngs/ngs_utils/ngsutils/bin/bamutils #version 0.5.9

module add picard-2.8.1
PICARD_RUN="java -Xmx12g -jar $PICARD281"
QUALIMAP=/storage/brno3-cerit/home/nikolato/ngs/ngs_utils/qualimap_v2.2.1/qualimap
MULTIQC=multiqc


SCORE=150 #FOR BWA
MAX_SOFTCLIP=0.20 # Maximal percentage of the reads allowed to be soft-clipped
MAX_PERC_OF_MM=0.10 # Max % of mismatches in the read

THREADS=9


#########################


### 1) extract the read ends - first and last 200 nucleotides 
cat $FASTQ/*.fastq > $ALL_READS
cutadapt -l 200 $ALL_READS > $FASTQ/merged_end.R1.fastq
cutadapt -l -200 $ALL_READS > $FASTQ/merged_end.R2.fastq


### 2) map reads on the masked reference
mkdir -p $ALIGNMENT_MASK
bwa mem -t $THREADS -T $SCORE -L 50 -v 1 -M $REF_DIR/CP004011.masked.fasta $FASTQ/merged_end.R1.fastq $FASTQ/merged_end.R2.fastq | samtools view -F 4 -Sb -h - > $ALIGNMENT_MASK/merged.tmp
samtools sort $ALIGNMENT_MASK/merged.tmp $ALIGNMENT_MASK/merged
samtools index $ALIGNMENT_MASK/merged.bam

samtools view -f 4 -Sb -h $ALIGNMENT_MASK/merged.tmp > $ALIGNMENT_MASK/merged.unmapped.bam
samtools index $ALIGNMENT_MASK/merged.unmapped.bam

samtools flagstat $ALIGNMENT_MASK/merged.bam > $ALIGNMENT_MASK/merged.bam.flagstat

samtools view -f66 $ALIGNMENT_MASK/merged.bam | cut -f 9 > $ALIGNMENT_MASK/merged.insert-sizes.txt

rm $ALIGNMENT_MASK/merged.tmp

### 3) READS FILTERING 
### 3_1) Too much softclipped
# Running it on the whole file and saving it creates error in BAM validation; we have to just take the read names from here and remove them from the alignment
$BAMUTILS removeclipping -f $ALIGNMENT_MASK/merged.bam $ALIGNMENT_MASK/merged.sfclip.bam.tmp
$SAMTOOLS view $ALIGNMENT_MASK/merged.sfclip.bam.tmp | grep 'ZC:f:' | awk '{for (i=1;i<=NF;i++){if ($i ~/ZC:f:/) {print $1, $i}}}' | sed 's/ZC:f://' | awk -v MAX_SOFTCLIP="$MAX_SOFTCLIP" ' $2 > MAX_SOFTCLIP {print $1}' > $ALIGNMENT_MASK/merged.bam.read_names_to_remove_highSoftClip.txt.tmp # Get only softclipped reads above the threshold
#rm $ALIGNMENT_MASK/merged.bam.tmp # Remove temporary BAMs
cat $ALIGNMENT_MASK/merged.bam.read_names_to_remove_highSoftClip.txt.tmp | cut -f 1 | sort -T $SCRATCH --parallel=$THREADS | uniq > $ALIGNMENT_MASK/merged.bam.read_names_to_remove_highSoftClip.txt # Get unique read names with softclipping!
#rm $ALIGNMENT_MASK/merged.bam.read_names_to_remove_highSoftClip.txt.tmp

# Use extracted filtered read name and filter it out from the original bam but ONLY when the soft clipping filtering had some results
if [ -s $ALIGNMENT_MASK/merged.bam.read_names_to_remove_highSoftClip.txt ]
then
	echo "Number of too much soft clipped reads is" `wc -l $ALIGNMENT_MASK/merged.bam.read_names_to_remove_highSoftClip.txt` "for sample " $ALIGNMENT_MASK/merged.bam.insert_size
	$PICARD_RUN FilterSamReads I=$ALIGNMENT_MASK/merged.bam O=$ALIGNMENT_MASK/merged.bam.filtered_out.bam READ_LIST_FILE=$ALIGNMENT_MASK/merged.bam.read_names_to_remove_highSoftClip.txt FILTER=includeReadList # Include reads
	$PICARD_RUN FilterSamReads I=$ALIGNMENT_MASK/merged.bam O=$ALIGNMENT_MASK/merged.sfclip.bam READ_LIST_FILE=$ALIGNMENT_MASK/merged.bam.read_names_to_remove_highSoftClip.txt FILTER=excludeReadList # Exclude reads
else
	echo "There are none to much soft clipped reads with " $MAX_SOFTCLIP " % soft clipping for sample " $ALIGNMENT_MASK/merged.bam. Continue without filtering.
	cp $ALIGNMENT_MASK/merged.bam $ALIGNMENT_MASK/merged.sfclip.bam
fi

$SAMTOOLS index $ALIGNMENT_MASK/merged.sfclip.bam

samtools flagstat $ALIGNMENT_MASK/merged.sfclip.bam > $ALIGNMENT_MASK/merged.bam.filter.flagstat
samtools view -f66 $ALIGNMENT_MASK/merged.sfclip.bam | cut -f 9 > $ALIGNMENT_MASK/merged_filter2.insert-sizes.txt	
echo "remove too much mismatches done"
#rm $ALIGNMENT_MASK/merged.bam.tmp

#############################
## 3_2) Remove too many mismatches
$BAMUTILS filter $ALIGNMENT_MASK/merged.sfclip.bam $ALIGNMENT_MASK/merged.sfclip.mm1.bam -failed $ALIGNMENT_MASK/merged.sfclip.mm.fail.txt -maximum_mismatch_ratio $MAX_PERC_OF_MM # Filter by perc. of mismatches; It's more filtering on edit distance = how many nucleotides have to be changed to get exactly the reference sequence; indels are counted as many times as they are "long"; error in the source https://github.com/ngsutils/ngsutils/pull/18
###$BAMUTILS filter $ALIGNMENT_MASK/merged.sfclip.mm1.bam $ALIGNMENT_MASK/merged.sfclip.mm2.bam -failed $ALIGNMENT_MASK/merged.sfclip.mm2.fail.txt -mismatch $MAX_NUMBER_OF_MM # Filter by num. of mismatches

###cat ${i%.*}.mm1.fail.txt ${i%.*}.mm2.fail.txt > ${i%.*}.mm.fail.txt # Merge mapping filtered by mismatches

# Use extracted filtered read name and filter it out from the original bam but ONLY when the mismatches filtering had some results
if [ -s $ALIGNMENT_MASK/merged.sfclip.mm.fail.txt ]
then
	cat $ALIGNMENT_MASK/merged.sfclip.mm.fail.txt | awk '{print $1}' | sort -T $SCRATCH --parallel=$THREADS | uniq > $ALIGNMENT_MASK/merged.sfclip.mm.fail.txt.tmp
	echo "Number of too much mismatched reads is" `wc -l $ALIGNMENT_MASK/merged.sfclip.mm.fail.txt.tmp` "for sample " $ALIGNMENT_MASK/merged.sfclip.bam
	$PICARD_RUN FilterSamReads I=$ALIGNMENT_MASK/merged.sfclip.bam O=$ALIGNMENT_MASK/merged.sfclip.mm.filtOut.bam READ_LIST_FILE=$ALIGNMENT_MASK/merged.sfclip.mm.fail.txt.tmp FILTER=includeReadList # Include reads
	#rm $ALIGNMENT_MASK/merged.sfclip.mm.fail.txt.tmp
else
	echo "There are none to much mismatched reads with " $MAX_PERC_OF_MM " % and " $MAX_NUMBER_OF_MM " mismatches for sample " $ALIGNMENT_MASK/merged.sfclip.bam". Nothing to report."
fi

samtools index $ALIGNMENT_MASK/merged.sfclip.mm1.bam
samtools flagstat $ALIGNMENT_MASK/merged.sfclip.mm1.bam > $ALIGNMENT_MASK/merged.sfclip.mm1.flagstat

#############################

###3_3) filter on RAW insert size
samtools view -h $ALIGNMENT_MASK/merged.sfclip.mm1.bam | \
awk 'substr($0,1,1)=="@" || ($9>=1000 && $9<=6000) || ($9<=-1000 && $9>=-6000)' > $ALIGNMENT_MASK/merged.sfclip.mm1.insert_size.tmp
samtools view -bS -h $ALIGNMENT_MASK/merged.sfclip.mm1.insert_size.tmp > $ALIGNMENT_MASK/merged.sfclip.mm1.insert_size.bam
samtools index $ALIGNMENT_MASK/merged.sfclip.mm1.insert_size.bam
samtools flagstat $ALIGNMENT_MASK/merged.sfclip.mm1.insert_size.bam > $ALIGNMENT_MASK/merged.sfclip.mm1.insert_size.flagstat
samtools view -f66 $ALIGNMENT_MASK/merged.sfclip.mm1.insert_size.bam | cut -f 9 > $ALIGNMENT_MASK/merged.sfclip.mm1.isert_size.insert-sizes.txt
echo "insert size filtering done"

###COVERAGE####
mosdepth $ALIGNMENT_MASK/merged.sfclip.mm1.insert_size $ALIGNMENT_MASK/merged.sfclip.mm1.insert_size.bam -b $BED
mosdepth $ALIGNMENT_MASK/merged.sfclip.mm1 $ALIGNMENT_MASK/merged.sfclip.mm1.bam -b $BED
mosdepth $ALIGNMENT_MASK/merged.sfclip $ALIGNMENT_MASK/merged.sfclip.bam -b $BED
mosdepth $ALIGNMENT_MASK/merged $ALIGNMENT_MASK/merged.bam -b $BED



### 4) SORT READS INTO GROUPS ACCORDING TO ROI IN BED FILE
mkdir -p $ALIGNMENT_MASK/bam_fastq
mkdir -p $LONG_READS

### create bam for each ROI
for bed_file in $TAG_DIR_MASKED/*.bed
do
	tag="$(basename $bed_file)"
	tag_name=${tag%*.bed}
	echo $tag_name
	samtools view -h -b $ALIGNMENT_MASK/merged.sfclip.mm1.insert_size.bam -L $bed_file > $ALIGNMENT_MASK/$tag_name.bam
	samtools sort $ALIGNMENT_MASK/$tag_name.bam $ALIGNMENT_MASK/$tag_name.sorted

	### get reads from bam file
	picard SamToFastq I=$ALIGNMENT_MASK/$tag_name.sorted.bam F=$ALIGNMENT_MASK/bam_fastq/$tag_name.R1.fastq F2=$ALIGNMENT_MASK/bam_fastq/$tag_name.R2.fastq

	### extract read names from from tags mapped to specific ROI and next extract long reads for each ROI
	grep -e ^@ $ALIGNMENT_MASK/bam_fastq/$tag_name.R1.fastq | grep /1$ | sed 's#/1$##g' | sed 's#^@##g' > $ALIGNMENT_MASK/bam_fastq/$tag_name.read_names
	filterbyname.sh in=$ALL_READS out=$LONG_READS/$tag_name.long_reads.fastq names=$ALIGNMENT_MASK/bam_fastq/$tag_name.read_names qin=33 overwrite=true include=t

done

### 5) FILTER READS specifically based on their expected length
cd $LONG_READS

for file in $LONG_READS/*.long_reads.fastq
do
	tag="$(basename $file)"
	tag_name=${tag%*.long_reads.fastq}

	echo $tag_name

	MIN=$(cat $INSERT_SIZE_FILE | grep $tag_name | awk '{print $3}')
	MAX=$(cat $INSERT_SIZE_FILE | grep $tag_name | awk '{print $4}')
	echo $MIN
	echo $MAX

	###keep only those long reads with appropriate read length
	cutadapt -m $MIN -M $MAX --too-short-output $tag_name.tooShort --too-long-output $tag_name.tooLong -o $tag_name.length_filter.fastq $file > $tag_name.cutadapt_report
done


### 6) DE NOVO ASSEMBLY

module add canu-1.8
mkdir -p $ASSEMBLY
mkdir -p $MAPPING
mkdir -p $RAW_REF_CONTIGS
mkdir -p $VARIANTS
mkdir -p $POLISHED_GENOME
mkdir -p $RACON
mkdir -p $MEDAKA
mkdir -p $TEST_RACON
mkdir -p $TEST_MEDAKA

for fastq in $LONG_READS/*.length_filter.fastq
do
	file="$(basename $fastq)"
	PREFIX=${file%*.length_filter.fastq}
	echo "processing tag $PREFIX"
	canu readSamplingCoverage=3000 maxThreads=8 -p $PREFIX -d $ASSEMBLY/$PREFIX.canu stopOnLowCoverage=1 genomeSize=$GENOME_SIZE useGrid=false -fast -nanopore-raw $fastq
	
 	###IF the longest read cover all the other reads, contigs are not created (file is empty) and will be represented by the longest read; If trimmed seqs are not generated, take the longest from corrected reads
	### canu, flye, unicycler, miniasm were tested, not working well. Problem with "short" long-read amplicons (they are uniform)

	if [ -s "$ASSEMBLY/$PREFIX.canu/$PREFIX.contigs.fasta" ]; 
	then
	echo "$ASSEMBLY/$PREFIX.canu/$PREFIX.contigs.fasta contain contigs and that's GOOOD" >> $RESULTS/contig_warnings.txt
	else
  	echo "$ASSEMBLY/$PREFIX.canu/$PREFIX.contigs.fasta is empty or missing. Contig will be represented by the longest trimmed sequence from $ASSEMBLY/$PREFIX.canu/$PREFIX.trimmedReads.fasta.gz!" >> $RESULTS/contig_warnings.txt
	echo "contigs from trimmed reads" > $ASSEMBLY/$PREFIX.canu/$PREFIX.WARNING.TXT
	gunzip -c $ASSEMBLY/$PREFIX.canu/$PREFIX.trimmedReads.fasta.gz | paste - - | awk -F '\t' '{L=length($2);if(L>M) {M=L;R=$0;}} END {print R;}' | tr "\t" "\n" > $ASSEMBLY/$PREFIX.canu/$PREFIX.contigs.fasta
	fi
	
	### if there are no contigs and trimmed reads neither, than use the longest corrected read as contigs
	if [ ! -s $ASSEMBLY/$PREFIX.canu/$PREFIX.trimmedReads.fasta.gz ]; 
	then
 	echo "$ASSEMBLY/$PREFIX.canu/$PREFIX.trimmedReads.fasta.gz is empty or missing. Contigs will be represented by the longest sequence from $ASSEMBLY/$PREFIX.canu/$PREFIX.correctedReads.fasta.gz!" >> $RESULTS/contig_warnings.txt
	echo "contigs from corrected reads" > $ASSEMBLY/$PREFIX.canu/$PREFIX.WARNING.TXT
	gunzip -c $ASSEMBLY/$PREFIX.canu/$PREFIX.correctedReads.fasta.gz | paste - - | awk -F '\t' '{L=length($2);if(L>M) {M=L;R=$0;}} END {print R;}' | tr "\t" "\n" > $ASSEMBLY/$PREFIX.canu/$PREFIX.contigs.fasta
	fi

	### if there are no corrected reads, contigs will not be present 
	if [ ! -s $ASSEMBLY/$PREFIX.canu/$PREFIX.correctedReads.fasta.gz ]; 
	then
 	echo "$ASSEMBLY/$PREFIX.canu/$PREFIX.correctedReads.fasta.gz is empty or missing. Contigs will NOT BE PRESENT!" >> $RESULTS/contig_warnings.txt
	echo "contigs for $ASSEMBLY/$PREFIX.canu/$PREFIX are NOT PRESENT" > $ASSEMBLY/$PREFIX.canu/$PREFIX.WARNING.TXT
	fi


###########################################################################################################################################################################################

	cp $ASSEMBLY/$PREFIX.canu/$PREFIX.contigs.fasta $RAW_REF_CONTIGS
	##rename contig header
	sed -i "1s/.*/\>$PREFIX.contigs.fasta/" $RAW_REF_CONTIGS/$PREFIX.contigs.fasta

	minimap2 -ax map-ont -t 4 $RAW_REF_CONTIGS/$PREFIX.contigs.fasta $fastq > $MAPPING/$PREFIX.sorted.sam
	samtools view -Sb -h $MAPPING/$PREFIX.sorted.sam > $MAPPING/$PREFIX.sorted.tmp
	samtools sort $MAPPING/$PREFIX.sorted.tmp $MAPPING/$PREFIX
	samtools index $MAPPING/$PREFIX.bam
	rm $MAPPING/*.tmp	
	
    ##COVERAGE
    mkdir -p $MAPPING/coverage
 	mosdepth $MAPPING/coverage/$PREFIX $MAPPING/$PREFIX.bam 
    echo "chrom	length	bases	mean	min	max" > $MAPPING/coverage/coverage_all_ROI.txt
    grep TP $MAPPING/coverage/*summary.txt >> $MAPPING/coverage/coverage_all_ROI.txt    ### doesn't contain genes/roi with zero coverage - see .cov for zero covered fragments
    bedtools genomecov -ibam $MAPPING/$PREFIX.bam -bga > $MAPPING/coverage/$PREFIX.cov

### 7) VARIANT CALLING

	###CALL VARIANTS WITH RACON

	
	racon -m 8 -x -6 -g -8 -w 500 -t 14 $fastq $MAPPING/$PREFIX.sorted.sam $RAW_REF_CONTIGS/$PREFIX.contigs.fasta > $RACON/$PREFIX.minimap.racon.fa
	minimap2 -ax map-ont -t 4 $RACON/$PREFIX.minimap.racon.fa $fastq > $RACON/$PREFIX.minimap_polish1.sam

	###ROUND 2-4
	racon -m 8 -x -6 -g -8 -w 500 -t 14 $fastq $RACON/$PREFIX.minimap_polish1.sam $RACON/$PREFIX.minimap.racon.fa > $RACON/$PREFIX.minimap.racon2.fa
	minimap2 -ax map-ont -t 4 $RACON/$PREFIX.minimap.racon2.fa $fastq > $RACON/$PREFIX.minimap_polish2.sam
	#--
	racon -m 8 -x -6 -g -8 -w 500 -t 14 $fastq $RACON/$PREFIX.minimap_polish2.sam $RACON/$PREFIX.minimap.racon2.fa > $RACON/$PREFIX.minimap.racon3.fa
	minimap2 -ax map-ont -t 4 $RACON/$PREFIX.minimap.racon3.fa $fastq > $RACON/$PREFIX.minimap_polish3.sam
	#--
	racon -m 8 -x -6 -g -8 -w 500 -t 14 $fastq $RACON/$PREFIX.minimap_polish3.sam $RACON/$PREFIX.minimap.racon3.fa > $RACON/$PREFIX.minimap.racon4.fa

	##CALL VARIANTS WITH MEDAKA
    conda activate medaka
    medaka_consensus -i $fastq -d $RAW_REF_CONTIGS/$PREFIX.contigs.fasta -o $MEDAKA/$PREFIX.medaka -t 4
    mv $MEDAKA/$PREFIX.medaka/consensus.fasta $MEDAKA/$PREFIX.medaka/$PREFIX.medaka.fasta

    medaka_consensus -i $fastq -d $RACON/$PREFIX.minimap.racon.fa -o $MEDAKA/$PREFIX.racon.medaka -t 4
    mv $MEDAKA/$PREFIX.racon.medaka/consensus.fasta $MEDAKA/$PREFIX.racon.medaka/$PREFIX.racon.medaka.fasta

    conda deactivate

    conda activate /storage/brno2/home/nikolato/.conda/envs/test


### 8) TEST PERFORMANCE
	samtools faidx $CONTROL_GENOME/$PREFIX.fa
	cp $RACON/$PREFIX*.racon*.fa $TEST_RACON
	cp -r $MEDAKA/$PREFIX.*medaka* $TEST_MEDAKA
	###TEST RACON
	cd $TEST_RACON
	echo "this is PREFIX $PREFIX"
    minimap2 -ax asm5 -t 4 $CONTROL_GENOME/$PREFIX.fa $PREFIX.minimap.racon4.fa > $PREFIX.minimap.racon4.sam 
    samtools view -Sb -h $PREFIX.minimap.racon4.sam > $PREFIX.minimap.racon4.sorted.tmp
    samtools sort $PREFIX.minimap.racon4.sorted.tmp $PREFIX.minimap.racon4
    samtools index $PREFIX.minimap.racon4.bam
    freebayes --bam $PREFIX.minimap.racon4.bam --fasta-reference $CONTROL_GENOME/$PREFIX.fa -C1 > $PREFIX.minimap.racon4.freebayes.vcf
    vcf2tsv $PREFIX.minimap.racon4.freebayes.vcf > $PREFIX.minimap.racon4.freebayes.tsv
    cat ./*tsv | grep TP > $SAMPLE.racon.allVars.tsv
	cd ..
	
	###TEST MEDAKA
	cd $TEST_MEDAKA
    minimap2 -ax asm5 -t 4 $CONTROL_GENOME/$PREFIX.fa $PREFIX.medaka/$PREFIX.medaka.fasta > $PREFIX.minimap.medaka.sam 
    samtools view -Sb -h $PREFIX.minimap.medaka.sam > $PREFIX.minimap.medaka.sorted.tmp
    samtools sort $PREFIX.minimap.medaka.sorted.tmp $PREFIX.minimap.medaka
    samtools index $PREFIX.minimap.medaka.bam
    freebayes --bam $PREFIX.minimap.medaka.bam --fasta-reference $CONTROL_GENOME/$PREFIX.fa -C1 > $PREFIX.minimap.medaka.freebayes.vcf
    vcf2tsv $PREFIX.minimap.medaka.freebayes.vcf > $PREFIX.minimap.medaka.freebayes.tsv
    cat ./*.minimap.medaka.freebayes.tsv | grep TP > $SAMPLE.medaka.allVars.tsv
    cd ..

    cd $TEST_MEDAKA
    minimap2 -ax asm5 -t 4 $CONTROL_GENOME/$PREFIX.fa $PREFIX.racon.medaka/$PREFIX.racon.medaka.fasta > $PREFIX.minimap.racon.medaka.sam 
    samtools view -Sb -h $PREFIX.minimap.racon.medaka.sam > $PREFIX.minimap.racon.medaka.sorted.tmp
    samtools sort $PREFIX.minimap.racon.medaka.sorted.tmp $PREFIX.minimap.racon.medaka
    samtools index $PREFIX.minimap.racon.medaka.bam
    freebayes --bam $PREFIX.minimap.racon.medaka.bam --fasta-reference $CONTROL_GENOME/$PREFIX.fa -C1 > $PREFIX.minimap.racon.medaka.freebayes.vcf
    vcf2tsv $PREFIX.minimap.racon.medaka.freebayes.vcf > $PREFIX.minimap.racon.medaka.freebayes.tsv
    cat ./*.minimap.racon.medaka.freebayes.tsv | grep TP > $SAMPLE.racon.medaka.allVars.tsv
    cd ..

done
