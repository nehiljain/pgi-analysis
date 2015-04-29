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

begin=$(date +%s)
$gatk -T CombineVariants -R $path_ref \
-V:pools $path_vcfC/pools_filtered_snps.vcf \
-V:indiv $path_vcfC/filtered_snps.vcf \
-genotypeMergeOptions UNIQUIFY \
-o $path_vcfC/cattle_combined_snps.vcf \
>$path_log/cattle.gatk.CombineVariants.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk CombineVariants runtime: $difftime seconds"

bgzip $path_vcfC/cattle_combined_snps.vcf
tabix $path_vcfC/cattle_combined_snps.vcf.gz

bcftools stats $path_vcfC/cattle_combined_snps.vcf.gz \
>$path_vcfC/cattle_combined.stats

gunzip $path_vcfC/cattle_combined_snps.vcf.gz

begin=$(date +%s)
$gatk -T SelectVariants -R $path_ref \
-V $path_vcfC/cattle_combined_snps.vcf \
-o $path_vcfC/cattle_biallelic_snps.vcf \
-restrictAllelesTo BIALLELIC \
>$path_log/cattle.gatk.SelectVariantsBiallelic.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk SelectVariants Biallelic runtime: $difftime seconds"

bgzip $path_vcfC/cattle_biallelic_snps.vcf
tabix $path_vcfC/cattle_biallelic_snps.vcf.gz

bcftools stats $path_vcfC/cattle_biallelic_snps.vcf.gz \
>$path_vcfC/cattle_biallelic.stats

gunzip $path_vcfC/cattle_biallelic_snps.vcf.gz

begin=$(date +%s)
$snpsift filter \
-f $path_vcfC/cattle_biallelic_snps.vcf \
"(QUAL >= 30) & ( AF < 1 ) & ( AF > 0 ) & ( DP >= 10 )" \
>$path_vcfC/cattle_filtered_snps.vcf \
2>$path_log/cattle.snpsift.snp.selection.log
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "snpsift snp filtration runtime: $difftime seconds" 

bgzip $path_vcfC/cattle_filtered_snps.vcf
tabix $path_vcfC/cattle_filtered_snps.vcf.gz

bcftools stats $path_vcfC/cattle_filtered_snps.vcf.gz \
>$path_vcfC/cattle_filtered.stats

gunzip $path_vcfC/cattle_filtered_snps.vcf.gz

echo ""
echo "------------------------------------------------------------"
