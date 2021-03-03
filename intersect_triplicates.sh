module add conda-modules-py37
conda activate /storage/brno2/home/nikolato/.conda/envs/test

RESULTS=/storage/brno3-cerit/home/nikolato/ngs/data/TPE_tags/OMPeome_8_10_2020/fastq/demultiplexed/results

ISEC_RACON=$RESULTS/isec_racon
ISEC_MEDAKA=$RESULTS/isec_medaka
ISEC_RACON_MEDAKA=$RESULTS/isec_racon_medaka

RACON_SUF=minimap.racon4.freebayes.norm.vcf.gz
MEDAKA_SUF=minimap.medaka.freebayes.norm.vcf.gz
RACON_MEDAKA_SUF=minimap.racon.medaka.freebayes.norm.vcf.gz


####################RUN INTERSECT FOR RACON
for i in {1,4,7}
do
    echo "triplicate ${i}"
    VARS_RACON=${RESULTS}/BC0${i}.pipeline4/racon/BC0${i}/performance_test
    VARS_RACON2=${RESULTS}/BC0$((i+1)).pipeline4/racon/BC0$((i+1))/performance_test
    VARS_RACON3=${RESULTS}/BC0$((i+2)).pipeline4/racon/BC0$((i+2))/performance_test
  
    echo "racon input folder: $VARS_RACON"
    ALL_SITES_RACON=$ISEC_RACON/BC0${i}.BC0$((i+1)).BC0$((i+2))/ALL_SITES_racon.txt

    for fastq in $VARS_RACON/*.minimap.racon4.freebayes.norm.vcf.gz
    do
        file="$(basename $fastq)"
        PREFIX=${file%*.minimap.racon4.freebayes.norm.vcf.gz}
        echo "processing tag $PREFIX"
        echo "processsing triplicate BC0${i}.BC0$((i+1)).BC0$((i+2))"

        echo "first file: $VARS_RACON/$PREFIX.$RACON_SUF"
        echo "second file: $VARS_RACON2/$PREFIX.$RACON_SUF"
        echo "third file: $VARS_RACON3/$PREFIX.$RACON_SUF"

        bcftools isec -p $ISEC_RACON/BC0${i}.BC0$((i+1)).BC0$((i+2))/$PREFIX -n +1 $VARS_RACON/$PREFIX.$RACON_SUF $VARS_RACON2/$PREFIX.$RACON_SUF $VARS_RACON3/$PREFIX.$RACON_SUF 
        
    done
   
    echo "results from intersecting are saved to: ${ALL_SITES_RACON}"
    cat $ISEC_RACON/BC0${i}.BC0$((i+1)).BC0$((i+2))/TP*/sites.txt > $ALL_SITES_RACON
done


#####################RUN INTERSECT FOR MEDAKA
for i in {1,4,7}
do
    echo "triplicate ${i}"
    VARS_MEDAKA=${RESULTS}/BC0${i}.pipeline4/medaka/BC0${i}/performance_test
    VARS_MEDAKA2=${RESULTS}/BC0$((i+1)).pipeline4/medaka/BC0$((i+1))/performance_test
    VARS_MEDAKA3=${RESULTS}/BC0$((i+2)).pipeline4/medaka/BC0$((i+2))/performance_test
  
    echo "medaka input folder: $VARS_MEDAKA"
    ALL_SITES_MEDAKA=$ISEC_MEDAKA/BC0${i}.BC0$((i+1)).BC0$((i+2))/ALL_SITES_medaka.txt

    for fastq in $VARS_MEDAKA/*.minimap.medaka.freebayes.norm.vcf.gz
    do
        file="$(basename $fastq)"
        PREFIX=${file%*.minimap.medaka.freebayes.norm.vcf.gz}
        echo "processing tag $PREFIX"
        echo "processsing triplicate BC0${i}.BC0$((i+1)).BC0$((i+2))"

        echo "first file: $VARS_MEDAKA/$PREFIX.$MEDAKA_SUF"
        echo "second file: $VARS_MEDAKA2/$PREFIX.$MEDAKA_SUF"
        echo "third file: $VARS_MEDAKA3/$PREFIX.$MEDAKA_SUF"

        bcftools isec -p $ISEC_MEDAKA/BC0${i}.BC0$((i+1)).BC0$((i+2))/$PREFIX -n +1 $VARS_MEDAKA/$PREFIX.$MEDAKA_SUF $VARS_MEDAKA2/$PREFIX.$MEDAKA_SUF $VARS_MEDAKA3/$PREFIX.$MEDAKA_SUF 
        
    done
   
    echo "results from intersecting are saved to: ${ALL_SITES_MEDAKA}"
    cat $ISEC_MEDAKA/BC0${i}.BC0$((i+1)).BC0$((i+2))/TP*/sites.txt > $ALL_SITES_MEDAKA
done

# #####################RUN INTERSECT FOR RACON_MEDAKA
for i in {1,4,7}
do
    echo "triplicate ${i}"
    VARS_RACON_MEDAKA=${RESULTS}/BC0${i}.pipeline4/medaka/BC0${i}/performance_test
    VARS_RACON_MEDAKA2=${RESULTS}/BC0$((i+1)).pipeline4/medaka/BC0$((i+1))/performance_test
    VARS_RACON_MEDAKA3=${RESULTS}/BC0$((i+2)).pipeline4/medaka/BC0$((i+2))/performance_test
  
    echo "racon.medaka input folder: $VARS_RACON_MEDAKA"
    ALL_SITES_RACON_MEDAKA=$ISEC_RACON_MEDAKA/BC0${i}.BC0$((i+1)).BC0$((i+2))/ALL_SITES_racon.medaka.txt

    for fastq in $VARS_RACON_MEDAKA/*.minimap.racon.medaka.freebayes.norm.vcf.gz
    do
        file="$(basename $fastq)"
        PREFIX=${file%*.minimap.racon.medaka.freebayes.norm.vcf.gz}
        echo "processing tag $PREFIX"
        echo "processsing triplicate BC0${i}.BC0$((i+1)).BC0$((i+2))"

        echo "first file: $VARS_RACON_MEDAKA/$PREFIX.$RACON_MEDAKA_SUF"
        echo "second file: $VARS_RACON_MEDAKA2/$PREFIX.$RACON_MEDAKA_SUF"
        echo "third file: $VARS_RACON_MEDAKA3/$PREFIX.$RACON_MEDAKA_SUF"

        bcftools isec -p $ISEC_RACON_MEDAKA/BC0${i}.BC0$((i+1)).BC0$((i+2))/$PREFIX -n +1 $VARS_RACON_MEDAKA/$PREFIX.$RACON_MEDAKA_SUF $VARS_RACON_MEDAKA2/$PREFIX.$RACON_MEDAKA_SUF $VARS_RACON_MEDAKA3/$PREFIX.$RACON_MEDAKA_SUF 
        
    done
   
    echo "results from intersecting are saved to: ${ALL_SITES_RACON_MEDAKA}"
    cat $ISEC_RACON_MEDAKA/BC0${i}.BC0$((i+1)).BC0$((i+2))/TP*/sites.txt > $ALL_SITES_RACON_MEDAKA
done
