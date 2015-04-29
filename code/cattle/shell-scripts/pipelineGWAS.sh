#$ -S /bin/bash
clear

## environment information
echo ""
echo ""
echo "============================================================"
echo "environment information:"
echo "Host: `hostname`"
echo "Login: $USER"
echo "Run time: `date`"
echo ""

#### loading cluster modules
module load picard/git
module load samtools/1.1

#### creating software enviroment
echo "============================================================"
echo "creating software environment:"
gatk="java -jar /home/kzukowski/soft/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar"
echo "gatk"
picard="java -jar /share/apps/picard/git/bin/picard.jar"
echo "picard"
echo "samtools"
snpeff="java -jar /home/kzukowski/soft/snpEff/snpEff.jar"
echo "snpeff"
snpsift="java -jar /home/kzukowski/soft/snpEff/SnpSift.jar"
echo "snpsift"
echo "bcftolls"
beagle="java -Xmx16g -jar /home/kzukowski/soft/beagle/beagle.r1399.jar"
echo "beagle"
beagle3="java -Xmx16g -jar /home/kzukowski/soft/beagle/beagle3.jar"
echo "beagle3"
plink="/home/kzukowski/soft/plink/plink-v1.903k"
echo "plink"
echo ""

#### creating paths
echo "============================================================"
echo "creating paths:"
echo ""
path_log="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_logs"; echo "path to log: $path_log"
path_vcf="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_ind_combined_vcf"; echo "path to vcf: $path_vcf"
path_csv="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_filtered_snp_list_csv"; echo "path to csv: $path_csv"
path_tmp="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/tmp"; echo "path to tmp: $path_tmp"
path_sync="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_pools_sync"; echo "path to sync: $path_sync"
path_syncF="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_filtered_sync"; echo "path to syncF: $path_syncF"
path_cmh="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_cmhtest_cmh"; echo "path to cmh: $path_cmh"
path_gwas="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_cmhtest_gwas"; echo "path to gwas: $path_gwas"
echo ""

echo "============================================================"
echo "analysis:"
echo "------------------------------------------------------------"

gawk '!/#/' $path_vcf/cattle_filtered_snps.vcf >$path_tmp/genome.step1
gawk '{print $1,$2}' $path_tmp/genome.step1 >$path_tmp/genome.step2

for i in `seq 1 29` X;
	do

	begin=$(date +%s)

	gawk '$1=="'$i'"' $path_tmp/genome.step2 >$path_tmp/$i.step3
	gawk '$1=$1' FS=" " OFS="," $path_tmp/$i.step3 >$path_tmp/$i.step4
	cat $path_csv/names $path_tmp/$i.step4 >$path_csv/snp_list_chr$i.csv

	rm -r $path_tmp/$i.step3
	rm -r $path_tmp/$i.step4

	gawk 'BEGIN { FS=","; OFS="\t" } {$1=$1; print}' $path_csv/snp_list_chr$i.csv >$path_tmp/snp_list_chr$i.tmp.tsv
	tail -n +2 $path_tmp/snp_list_chr$i.tmp.tsv >$path_tmp/snp_list_chr$i.tsv

	gawk '{OFS="\t"} FNR==NR{a[$2]=$2;next}(a[$2]){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19, a[$1]}' \
	$path_tmp/snp_list_chr$i.tsv $path_sync/$i.sync >$path_tmp/chr$i.filtered.tmp.sync

	echo $i
	echo "snps in raw sync "; wc  -l $path_sync/$i.sync
	echo "snps in list: "; wc -l $path_tmp/snp_list_chr$i.tsv
	echo "snps in filtered sync: "; wc -l $path_tmp/chr$i.filtered.tmp.sync

	sed 's/\t$//' $path_tmp/chr$i.filtered.tmp.sync >$path_syncF/chr$i.filtered.sync

	rm -r $path_tmp/snp_list_chr$i.tmp.tsv
	rm -r $path_tmp/snp_list_chr$i.tsv
	rm -r $path_tmp/chr$i.filtered.tmp.sync

	perl /home/kzukowski/soft/popoolation2_1201/cmh-test.pl \
	--input $path_syncF/chr$i.filtered.sync \
	--output $path_cmh/chr$i.cmh \
	--min-count 2 \
	--min-coverage 2 \
	--max-coverage 1000000000 \
	--population 9-1,10-2,11-3,12-4,13-5,14-6,15-7,16-8 \
	--remove-temp \
	>$path_log/popoolation2.chr$i.log 2>&1

	echo "snps in cmh: "; wc -l $path_cmh/chr$i.cmh

	cut -f 1,2,3,20 $path_cmh/chr$i.cmh >$path_gwas/chr$i.gwas

	echo "snps in gwas: "; wc -l $path_gwas/chr$i.gwas

	termin=$(date +%s)
	difftime=$(( $termin - $begin ))
	echo "popoolation2 runtime: $difftime seconds"
	echo "------------------------------------------------------------"

	done

rm -r $path_tmp/genome.step1
rm -r $path_tmp/genome.step2

echo ""
echo "============================================================"
