
#!/bin/bash
# set -euo pipefail
#PBS -l walltime=40:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

# Description: Get FASTQ metrics for validation run .
# Authors: Christopher Medway and Laura McCluskey, Seemu Ali
# Date: 27th March 2020
# Usage: qsub run_star-fusion.sh [inside sample dir with .variables and \
# .fastq.gz files]

version=0.0.4

# source variables file
. *.variables

mkdir /data/results/$seqId/$panel/Summary_Metrics
# copy the panel & pipeline variables locally and source
cp /data/diagnostics/pipelines/SomaticFusionHaem/SomaticFusionHaem-$version/$panel/$panel.variables . 
cp /data/diagnostics/pipelines/SomaticFusionHaem/SomaticFusionHaem-$version/SomaticFusion.config .
cp /data/diagnostics/pipelines/SomaticFusionHaem/SomaticFusionHaem-$version/make-fusion-report.py .
cp /data/diagnostics/pipelines/SomaticFusionHaem/SomaticFusionHaem-$version/RNA_fusion_group_file.txt .
cp /data/diagnostics/pipelines/SomaticFusionHaem/SomaticFusionHaem-$version/RNAFusion-ROI_adapted.bed .
cp /data/diagnostics/pipelines/SomaticFusionHaem/SomaticFusionHaem-$version/fusion_report_referrals.py .
cp /data/diagnostics/pipelines/SomaticFusionHaem/SomaticFusionHaem-$version/fusion_report_cases.py . 
cp /data/diagnostics/pipelines/SomaticFusionHaem/SomaticFusionHaem-$version/total_reads_list.py /data/results/$seqId/$panel/Summary_Metrics
gatk3=/share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar
minMQS=20
minBQS=10

vendorCaptureBed=./180702_HG19_PanCancer_EZ_capture_targets.bed
vendorPrimaryBed=./180702_HG19_PanCancer_EZ_primary_targets.bed

. $panel.variables
. SomaticFusion.config

# set conda env
source "$conda_bin_path"/activate SomaticFusion


#######################
# Preprocessing FASTQ #
#######################

#count how many core FASTQC tests failed
countQCFlagFails() {
    grep -E "Basic Statistics|Per base sequence quality|Per tile sequence quality|Per sequence quality scores|Per base N content" "$1" | \
    grep -v ^PASS | \
    grep -v ^WARN | \
    wc -l | \
    sed 's/^[[:space:]]*//g'
}

#record FASTQC pass/fail
rawSequenceQuality=PASS



for fastqPair in $(ls "$sampleId"_S*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do

    #parse fastq filenames
    laneId=$(echo "$fastqPair" | cut -d_ -f3)
    read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
    read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)

    #trim adapters
    cutadapt -a $read1Adapter -A $read2Adapter -m 35 -o "$seqId"_"$sampleId"_"$laneId"_R1.fastq -p "$seqId"_"$sampleId"_"$laneId"_R2.fastq "$read1Fastq" "$read2Fastq"


    fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R1.fastq
    fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R2.fastq

    mv "$seqId"_"$sampleId"_"$laneId"_R1_fastqc/summary.txt "$seqId"_"$sampleId"_"$laneId"_R1_fastqc.txt
    mv "$seqId"_"$sampleId"_"$laneId"_R2_fastqc/summary.txt "$seqId"_"$sampleId"_"$laneId"_R2_fastqc.txt

    #check FASTQC output
    if [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R1_fastqc.txt) -gt 0 ] || [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R2_fastqc.txt) -gt 0 ]; then
        rawSequenceQuality=FAIL
    fi


############
# Run STAR #
############


STAR --chimSegmentMin 12 \
     --chimJunctionOverhangMin 12 \
     --chimSegmentReadGapMax 3 \
     --alignSJDBoverhangMin 10 \
     --alignMatesGapMax 200000 \
     --chimOutType WithinBAM \
     --alignIntronMax 200000  \
     --alignSJstitchMismatchNmax 5 -1 5 5  \
     --twopassMode Basic \
     --outSAMtype BAM Unsorted \
     --readFilesIn "$seqId"_"$sampleId"_"$laneId"_R1.fastq "$seqId"_"$sampleId"_"$laneId"_R2.fastq \
     --genomeDir /share/apps/star-fusion/GRCh37_gencode_v19_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx  \
     --outFileNamePrefix "$seqId"_"$sampleId"_"$laneId"_ \
     --outSAMattrRGline ID:"$sampleId" SM:"$sampleId"
done


###################
# Run STAR-Fusion #
###################


## STAR-fusion required a samples_file in order to handle multi-lane runs
R1=$(for i in ./"$sampleId"_*.fastq.gz; do echo $i | grep "R1"; done)
R2=$(for i in ./"$sampleId"_*.fastq.gz; do echo $i | grep "R2"; done)
SAMPLE=$(for i in $R1;do echo $sampleId;done)
paste <(printf %s "$SAMPLE") <(printf %s "$R1") <(printf %s "$R2") > "$sampleId".samples

# run STAR-Fusion
STAR-Fusion --genome_lib_dir $starfusion_lib \
            --samples_file $sampleId.samples \
            --output_dir ./STAR-Fusion/ \
            --FusionInspector validate \
            --denovo_reconstruct \
            --examine_coding_effect \
            --CPU $ncpus \
            --min_FFPM 1



############################
# merge and sort bam files #
############################

#samtools merge "$sampleId"_Aligned_out.bam *Aligned.out.bam #ONE LANE ONLY 

samtools sort "$seqId"_"$sampleId"_"$laneId"_Aligned.out.bam -o "$sampleId"_Aligned_sorted.bam

samtools index "$sampleId"_Aligned_sorted.bam

################################
# Generate Star Fusion Reports #
################################



python make-fusion-report.py \
    --sampleId $sampleId \
    --seqId $seqId \
    --panel $panel \
    --ip $(hostname --ip-address)

#deactivate conda env
source "$conda_bin_path"/deactivate



##############
# Run arriba #
##############
mkdir Arriba 

source "$conda_bin_path"/activate SomaticFusion

arriba -x "$sampleId"_Aligned_sorted.bam \
-g /share/apps/star-fusion/GRCh37_gencode_v19_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf \
-a /share/apps/star-fusion/GRCh37_gencode_v19_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa \
-o ./Arriba/"$sampleId"_arriba_fusions.tsv \
-b /home/transfer/miniconda3/envs/arriba/var/lib/arriba/blacklist_hg19_hs37d5_GRCh37_2018-11-04.tsv.gz \
-O ./Arriba/"$sampleId"_fusions_arriba_discarded.tsv \
-R 0 \
-f intragenic_exonic,same_gene

source "$conda_bin_path"/deactivate



#######################################
#Separate reports into referral types #
#######################################


mkdir ./Results
mkdir ./Results/Full_Reports 

source "$conda_bin_path"/activate SomaticFusion

python fusion_report_cases.py /data/results/"$seqId"/"$panel"/"$sampleId" /data/diagnostics/pipelines/SomaticFusionHaem/SomaticFusionHaem-$version/$panel $sampleId

source "$conda_bin_path"/deactivate



###############################
# Calculate depth of coverage #
###############################


/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar $gatk3 \
   -T DepthOfCoverage \
   -R /share/apps/star-fusion/star-fusion-ref.fa \
   -I "$sampleId"_Aligned_sorted.bam \
   -L RNAFusion-ROI_adapted.bed \
   -o "$seqId"_"$sampleId"_DepthOfCoverage \
   --countType COUNT_FRAGMENTS \
   --minMappingQuality $minMQS \
   --minBaseQuality $minBQS \
   -ct 30  \
   --omitLocusTable \
   -dt NONE \
   -U ALLOW_N_CIGAR_READS \
   -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
   -rf MappingQualityUnavailable



############################
# Run CoverageCalculatorPy #
############################

sed 's/:/\t/g'  "$seqId"_"$sampleId"_DepthOfCoverage | grep -v 'Locus' | sort -k1,1 -k2,2n | bgzip > "$seqId"_"$sampleId"_DepthOfCoverage.gz
tabix -b 2 -e 2 -s 1 "$seqId"_"$sampleId"_DepthOfCoverage.gz


source "$conda_bin_path"/activate CoverageCalculatorPy


python /home/transfer/pipelines/CoverageCalculatorPy/CoverageCalculatorPy.py \
-B RNAFusion-ROI_adapted.bed \
-D /data/results/"$seqId"/"$panel"/"$sampleId"/"$seqId"_"$sampleId"_DepthOfCoverage.gz \
--padding 0 \
--groupfile RNA_fusion_group_file.txt \
--outname "$sampleId"_coverage 



########################
# Calculate qc metrics #
########################


if [ ! -f ../"$seqId"_"$sampleId"_HsMetrics.txt ]; then

/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g \
     -jar /share/apps/picard-tools-distros/picard-tools-2.18.5/picard.jar CollectHsMetrics \
     I=/data/results/"$seqId"/"$panel"/"$sampleId"/"$sampleId"_Aligned_sorted.bam \
     O=/data/diagnostics/pipelines/SomaticFusionHaem/SomaticFusionHaem-$version/$panel/Summary_Metrics//"$seqId"_"$sampleId"_HsMetrics.txt \
     R=/share/apps/star-fusion/GRCh37_gencode_v19_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa \
     BAIT_INTERVALS=/data/diagnostics/pipelines/SomaticFusionHaem/SomaticFusionHaem-$version/$panel/trusight-rna-fusion.targeted_regions.sorted.bed \
     TARGET_INTERVALS=/data/diagnostics/pipelines/SomaticFusionHaem/SomaticFusionHaem-$version/$panel/trusight-rna-fusion.targeted_regions.sorted.bed \
     MAX_RECORDS_IN_RAM=2000000 \
     TMP_DIR=/state/partition1/tmpdir \
     MINIMUM_MAPPING_QUALITY=$minMQS \
     MINIMUM_BASE_QUALITY=$minBQS \
     CLIP_OVERLAPPING_READS=false

fi



#Alignment metrics: library sequence similarity
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g \
    -jar /share/apps/picard-tools-distros/picard-tools-2.18.5/picard.jar CollectAlignmentSummaryMetrics \
    R=/share/apps/star-fusion/GRCh37_gencode_v19_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa \
    I="$sampleId"_Aligned_sorted.bam \
    O=/data/results/$seqId/$panel/Summary_Metrics/"$sampleId"_AlignmentSummaryMetrics.txt \
    MAX_RECORDS_IN_RAM=2000000 \
TMP_DIR=/state/partition1/tmpdir


source "$conda_bin_path"/deactivate




#####################
# remove duplicates #
#####################

#remove duplicates

/share/apps/jre-distros/jre1.8.0_131/bin/java \
    -XX:GCTimeLimit=50 \
    -XX:GCHeapFreeLimit=10 \
    -Djava.io.tmpdir=/state/partition1/tmpdir \
    -Xmx2g \
    -jar /share/apps/picard-tools-distros/picard-tools-2.18.5/picard.jar \
    MarkDuplicates \
    I="$sampleId"_Aligned_sorted.bam \
    OUTPUT="$sampleId"_rmdup.bam \
    METRICS_FILE="$sampleId"_markDuplicatesMetrics.txt \
    CREATE_INDEX=true \
    MAX_RECORDS_IN_RAM=2000000 \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=/state/partition1/tmpdir \
    QUIET=true \
    VERBOSITY=ERROR \
    REMOVE_DUPLICATES=TRUE



##########################
# recalculate qc metrics #
##########################



#Alignment metrics: library sequence similarity
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g \
    -jar /share/apps/picard-tools-distros/picard-tools-2.18.5/picard.jar CollectAlignmentSummaryMetrics \
    R=/share/apps/star-fusion/GRCh37_gencode_v19_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa \
    I="$sampleId"_rmdup.bam \
    O=/data/results/$seqId/$panel/Summary_Metrics/"$sampleId"_AlignmentSummaryMetrics_rmdup.txt \
    MAX_RECORDS_IN_RAM=2000000 \
TMP_DIR=/state/partition1/tmpdir


##########################
# recalculate Coverage   #
##########################


#Calculate depth of coverage

/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar $gatk3 \
   -T DepthOfCoverage \
   -R /share/apps/star-fusion/star-fusion-ref.fa \
   -I "$sampleId"_rmdup.bam \
   -L RNAFusion-ROI_adapted.bed \
   -o "$sampleId"_rmdup_DepthOfCoverage \
   --countType COUNT_FRAGMENTS \
   --minMappingQuality 20  \
   --minBaseQuality 10 \
   -ct 30  \
   --omitLocusTable \
   -dt NONE \
   -U ALLOW_N_CIGAR_READS \
   -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
   -rf MappingQualityUnavailable



#run coverageCalculatorPy

source /home/transfer/miniconda3/bin/activate CoverageCalculatorPy

sed 's/:/\t/g'  "$sampleId"_rmdup_DepthOfCoverage | grep -v 'Locus' | sort -k1,1 -k2,2n | bgzip > "$sampleId"_rmdup_DepthOfCoverage.gz

tabix -b 2 -e 2 -s 1 "$sampleId"_rmdup_DepthOfCoverage.gz



python /home/transfer/pipelines/CoverageCalculatorPy/CoverageCalculatorPy.py \
-B RNAFusion-ROI_adapted.bed \
-D /data/results/"$seqId"/"$panel"/"$sampleId"/"$sampleId"_rmdup_DepthOfCoverage.gz \
--padding 0 \
--groupfile RNA_fusion_group_file.txt \
--outname "$sampleId"_rmdup_coverage


source "$conda_bin_path"/deactivate




###########################
# Generate Metrics Report #
###########################



source "$conda_bin_path"/activate SomaticFusion


if [ -e /data/results/$seqId/$panel/Summary_Metrics/"$sampleId"_AlignmentSummaryMetrics_rmdup.txt ]
then
    echo $sampleId >> /data/results/$seqId/$panel/Summary_Metrics/samples_list.txt
fi

numberSamplescomplete=$(cat /data/results/$seqId/$panel/Summary_Metrics/samples_list.txt | uniq | wc -l)
numberSamplesInProject=$(find ../ -maxdepth 2 -mindepth 2 | grep .variables | grep -v "$panel".variables | uniq | wc -l)

echo $numberSamplescomplete
echo $numberSamplesInProject 
find ../ -maxdepth 2 -mindepth 2 | grep .variables | grep -v "$panel".variables > proj_samp.txt
cat /data/results/$seqId/$panel/Summary_Metrics/samples_list.txt | uniq > smpcomp.txt




if [ $numberSamplescomplete -eq $numberSamplesInProject ]
then

    python /data/results/$seqId/$panel/Summary_Metrics/total_reads_list.py "$seqId" "/data/results/$seqId/$panel/Summary_Metrics"

else
    echo "not all samples have completed running. Finishing process for this sample."
fi



source "$conda_bin_path"/deactivate