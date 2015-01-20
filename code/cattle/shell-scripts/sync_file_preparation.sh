#!/bin/bash
clear



echo ""
echo "Host: `hostname`"
echo "Login: $USER"
echo "Run time: `date`"
echo ""

date
echo
echo "picard CleanSam"
echo

java -Xmx8g -jar /home/kzukowski/soft/picard-tools-1.119/CleanSam.jar \
INPUT=/share/volatile_scratch/nehil/cattle/Bos_taurus_DNAalignment_975/12766/run2332_4/12766.MPS12301871-B04.sorted.bam \
OUTPUT=/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.CleanSam.bam \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true \
TMP_DIR=/tmp \
>/share/volatile_scratch/kzukowski/pgi/cattle/test/log/12766.picard.CleanSam.log 2>&1

date
echo
echo "samtools view"
echo

samtools view -bq 20 \
/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.CleanSam.bam \
>/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.MAPQ20.bam \
2>/share/volatile_scratch/kzukowski/pgi/cattle/test/log/12766.samtools.MAPQ20.log

date
echo
echo "picard FixMateInformation"
echo

java -Xmx8g -jar /home/kzukowski/soft/picard-tools-1.119/FixMateInformation.jar \
INPUT=/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.MAPQ20.bam \
OUTPUT=/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.FixMate.bam \
SORT_ORDER=coordinate \
ASSUME_SORTED=true \
ADD_MATE_CIGAR=true \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true \
TMP_DIR=/tmp \
>/share/volatile_scratch/kzukowski/pgi/cattle/test/log/12766.picard.FixMate.log 2>&1

date
echo
echo "picard MarkDuplicates"
echo

java -Xmx8g -jar /home/kzukowski/soft/picard-tools-1.119/MarkDuplicates.jar \
INPUT=/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.FixMate.bam \
OUTPUT=/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.dedup.bam \
METRICS_FILE=/share/volatile_scratch/kzukowski/pgi/cattle/test/stats/12766.metrics.file \
REMOVE_DUPLICATES=true \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true \
TMP_DIR=/tmp \
>/share/volatile_scratch/kzukowski/pgi/cattle/test/log/12766.picard.MarkDuplicates_aln.log 2>&1

date
echo
echo "picard CollectMultipleMetrics"
echo

java -Xmx8g -jar /home/kzukowski/soft/picard-tools-1.119/CollectMultipleMetrics.jar \
INPUT=/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.dedup.bam \
OUTPUT=/share/volatile_scratch/kzukowski/pgi/cattle/test/stats/12766.CollectMultipleMetrics \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=/tmp \
>/share/volatile_scratch/kzukowski/pgi/cattle/test/log/12766.picard.CollectMultipleMetrics.log 2>&1

date
echo
echo "picard CollectGcBiasMetrics"
echo

java -Xmx8g -jar /home/kzukowski/soft/picard-tools-1.119/CollectGcBiasMetrics.jar \
INPUT=/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.dedup.bam \
OUTPUT=/share/volatile_scratch/kzukowski/pgi/cattle/test/stats/12766.CollectGcBiasMetrics \
CHART_OUTPUT=/share/volatile_scratch/kzukowski/pgi/cattle/test/stats/12766.CollectGcBiasMetrics.pdf \
REFERENCE_SEQUENCE=/share/volatile_scratch/kzukowski/pgi/cattle/reference/Bos_taurus.UMD3.1.dna.toplevel.fa \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=/tmp \
>/share/volatile_scratch/kzukowski/pgi/cattle/test/log/12766.picard.CollectGcBiasMetricss.log 2>&1

date
echo
echo "samtools mpileup"
echo



# merge all the bam files to mpileup
#
## -r 1 -> only for chr1

# for 1:X,Y

samtools mpileup -B -Q 20 -q 20 -r 1 \
-f /share/volatile_scratch/kzukowski/pgi/cattle/reference/Bos_taurus.UMD3.1.dna.toplevel.fa \
/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.dedup.bam \
>/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.mpileup \
2>/share/volatile_scratch/kzukowski/pgi/cattle/test/log/12766.samtools.mpileup.log





date
echo
echo "popoolation2/mpileup2sync"
echo

java -Xmx8g -jar /home/kzukowski/soft/popoolation2_1201/mpileup2sync.jar \
--input /share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.mpileup \
#--input next... and so on\
--output /share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.sync \
--fastq-type sanger \
--min-qual 20 \
--threads 2 \
>/share/volatile_scratch/kzukowski/pgi/cattle/test/log/12766.mpileup2sync.log 2>&1



