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
path_ref="/share/volatile_scratch/nehil/pgi_wc/cattle/reference/77/Bos_taurus.UMD3.1.dna.toplevel.fa"; echo "path to reference genome: $path_ref"
path_log="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_logs"; echo "path to log: $path_log"
path_tmp="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/tmp"; echo "path to tmp: $path_tmp"
path_stats="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_ind_stats"; echo "path to stats: $path_stats"
path_vcf0="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/gq_pools_snpcalls_freebayes_vcf"; echo "path to vcf0: $path_vcf0"
path_vcf="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_ind_snpcalls_sorted_vcf"; echo "path to vcf: $path_vcf"
path_vcfC="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_ind_combined_vcf"; echo "path to combined vcf: $path_vcfC"
echo ""

echo "============================================================"
echo "analysis:"


for i in PGI_Pool17 PGI_Pool19 PGI_Pool21 PGI_Pool23 PGI_Pool25 PGI_Pool27 PGI_Pool29 PGI_Pool31 PGI_Pool18 PGI_Pool20 PGI_Pool22 PGI_Pool24 PGI_Pool26 PGI_Pool28 PGI_Pool30 PGI_Pool32;
	do

	bgzip $path_vcf0/$i.vcf
	tabix $path_vcf0/$i.vcf.gz

	bcftools stats $path_vcf0/$i.vcf.gz >$path_vcf0/$i.stats

	gunzip $path_vcf0/$i.vcf.gz

	begin=$(date +%s)
	$picard SortVcf \
	I=$path_vcf0/$i.vcf \
	O=$path_vcf/$i.sort.vcf \
	SEQUENCE_DICTIONARY=/share/volatile_scratch/nehil/pgi_wc/cattle/reference/77/Bos_taurus.UMD3.1.dna.toplevel.dict \
	>$path_log/$i.picard.SortVcf.log 2>&1
	termin=$(date +%s)
	difftime=$(( $termin - $begin ))
	echo "picard SortVcf sample $i runtime: $difftime seconds"

	done

begin=$(date +%s)
$gatk -T CombineVariants -R $path_ref \
-V:Pool25 $path_vcf/PGI_Pool25.sort.vcf \
-V:Pool26 $path_vcf/PGI_Pool26.sort.vcf \
-V:Pool27 $path_vcf/PGI_Pool27.sort.vcf \
-V:Pool28 $path_vcf/PGI_Pool28.sort.vcf \
-V:Pool29 $path_vcf/PGI_Pool29.sort.vcf \
-V:Pool30 $path_vcf/PGI_Pool30.sort.vcf \
-V:Pool31 $path_vcf/PGI_Pool31.sort.vcf \
-V:Pool32 $path_vcf/PGI_Pool32.sort.vcf \
-V:Pool17 $path_vcf/PGI_Pool17.sort.vcf \
-V:Pool18 $path_vcf/PGI_Pool18.sort.vcf \
-V:Pool19 $path_vcf/PGI_Pool19.sort.vcf \
-V:Pool20 $path_vcf/PGI_Pool20.sort.vcf \
-V:Pool21 $path_vcf/PGI_Pool21.sort.vcf \
-V:Pool22 $path_vcf/PGI_Pool22.sort.vcf \
-V:Pool23 $path_vcf/PGI_Pool23.sort.vcf \
-V:Pool24 $path_vcf/PGI_Pool24.sort.vcf \
-genotypeMergeOptions UNIQUIFY \
-o $path_vcfC/pools_combined.vcf \
>$path_log/pools_gatk.CombineVariants.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk CombineVariants runtime: $difftime seconds"

begin=$(date +%s)
$gatk -T SelectVariants -R $path_ref \
-V $path_vcfC/pools_combined.vcf \
-selectType SNP \
-o $path_vcfC/pools_raw_snps.vcf \
>$path_log/pools_gatk.SelectVariantsSNP.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk SelectVariants SNP runtime: $difftime seconds"

begin=$(date +%s)
$snpsift filter \
-f $path_vcfC/pools_raw_snps.vcf \
"(QUAL >= 30) & ( AF < 1 ) & ( AF > 0 ) & ( DP >= 10 )" \
>$path_vcfC/pools_filtered_snps.vcf \
2>$path_log/pools.snpsift.snp.selection.log
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "snpsift snp filtration runtime: $difftime seconds" 

begin=$(date +%s)
$gatk -T SelectVariants -R $path_ref \
-V $path_vcfC/pools_combined.vcf \
-selectType INDEL \
-o $path_vcfC/pools_raw_indels.vcf \
>$path_log/pools_gatk.SelectVariantsINDEL.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk SelectVariants INDEL runtime: $difftime seconds"

begin=$(date +%s)
$snpsift filter \
-f $path_vcfC/pools_raw_indels.vcf \
"(QUAL >= 30) & ( AF < 1 ) & ( AF > 0 ) & ( DP >= 10 )" \
>$path_vcfC/pools_filtered_indels.vcf \
2>$path_log/pools.snpsift.indel.selection.log
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "snpsift indel filtration runtime: $difftime seconds"

echo ""
echo "------------------------------------------------------------"
