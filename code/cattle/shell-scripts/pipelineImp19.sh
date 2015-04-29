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
path_vcf0="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/gq_snpcalls_vcf"; echo "path to raw vcf: $path_vcf0"
path_vcf="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_ind_snpcalls_sorted_vcf"; echo "path to vcf: $path_vcf"
path_vcfC="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_ind_combined_vcf"; echo "path to combined vcf: $path_vcfC"
path_imp="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_ind_imputation"; echo "path to imputation vcf: $path_imp"
path_cc="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_ind_casecontol_vcf"; echo "path to cc vcf: $path_cc"
path_ld="/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/kacper_ind_ld"; echo "path to linkage disequilibrium vcf: $path_ld"
echo ""

echo "============================================================"
echo "analysis:"

for i in 19;
	do

	mkdir $path_imp/chr$i 

	begin=$(date +%s)
	$gatk -T SelectVariants -R $path_ref \
	-V $path_vcfC/filtered_snps.vcf \
	-ef \
	-L $i \
	-o $path_imp/chr$i/chr$i.filtered_snps.PASS.vcf \
	>$path_log/chr$i.gatk.SelectVariants.log 2>&1
	termin=$(date +%s)
	difftime=$(( $termin - $begin ))
	echo "gatk SelectVariants chromosome $i runtime: $difftime seconds"

	begin=$(date +%s)
	$snpsift filter \
	-f $path_imp/chr$i/chr$i.filtered_snps.PASS.vcf \
	"( AN >= 50 ) & ( AC >= 50 ) & ( AF < 1 ) & ( AF > 0 ) & ( DP >= 10 )" \
	>$path_imp/chr$i/chr$i.filtered_snps.selected.vcf \
	2>$path_log/chr$i.snpsift.selection.log
	termin=$(date +%s)
	difftime=$(( $termin - $begin ))
	echo "snpsift filtration runtime: $difftime seconds"

	bgzip $path_imp/chr$i/chr$i.filtered_snps.selected.vcf
	tabix $path_imp/chr$i/chr$i.filtered_snps.selected.vcf.gz

	bcftools stats $path_imp/chr$i/chr$i.filtered_snps.selected.vcf.gz \
	>$path_imp/chr$i/chr$i.stats

	begin=$(date +%s)
	$beagle \
	gt=$path_imp/chr$i/chr$i.filtered_snps.selected.vcf.gz \
	ped=/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/pedigree/pedigree \
	burnin-its=100 phase-its=1000 impute-its=1000 \
	out=$path_imp/chr$i/chr$i.imputed \
	>$path_log/chr$i.beagle.log 2>&1
	termin=$(date +%s)
	difftime=$(( $termin - $begin ))
	echo "beagle imputation runtime: $difftime seconds"

	gunzip $path_imp/chr$i/chr$i.imputed.vcf.gz

	bgzip $path_imp/chr$i/chr$i.imputed.vcf
	tabix $path_imp/chr$i/chr$i.imputed.vcf.gz

	bcftools stats $path_imp/chr$i/chr$i.imputed.vcf.gz \
	>$path_imp/chr$i/chr$i.stats.imputed

	gunzip $path_imp/chr$i/chr$i.filtered_snps.selected.vcf.gz
	gunzip $path_imp/chr$i/chr$i.imputed.vcf.gz

	begin=$(date +%s)
	$snpsift caseControl \
	-tfam $path_cc/samples.tped \
	$path_imp/chr$i/chr$i.imputed.vcf \
	>$path_cc/chr$i.cc.vcf \
	2>$path_log/chr$i.snpsift.caseControl.log
	termin=$(date +%s)
	difftime=$(( $termin - $begin ))
	echo "snpsift caseControl runtime: $difftime seconds"

	gunzip $path_imp/chr$i/chr$i.filtered_snps.selected.vcf.gz

	mkdir $path_ld/chr$i
	mkdir $path_ld/chr$i\i

	$plink --cow --vcf $path_imp/chr$i/chr$i.filtered_snps.selected.vcf --r2 --ld-window-r2 0 --ld-window-kb 1000 --ld-window 1000 --out $path_ld/chr$i/ld
	$plink --cow --vcf $path_imp/chr$i/chr$i.filtered_snps.selected.vcf --blocks no-pheno-req no-small-max-span --blocks-max-kb 1000 --out $path_ld/chr$i/blocks
	$plink --cow --vcf $path_imp/chr$i/chr$i.filtered_snps.selected.vcf --show-tags all --tag-kb 1000 --out $path_ld/chr$i/tags

	$plink --cow --vcf $path_imp/chr$i/chr$i.imputed.vcf --r2 --ld-window-r2 0 --ld-window-kb 1000 --ld-window 1000 --out $path_ld/chr$i\i/ld
	$plink --cow --vcf $path_imp/chr$i/chr$i.imputed.vcf --blocks no-pheno-req no-small-max-span --blocks-max-kb 1000 --out $path_ld/chr$i\i/blocks
	$plink --cow --vcf $path_imp/chr$i/chr$i.imputed.vcf --show-tags all --tag-kb 1000 --out $path_ld/chr$i\i/tags

	done

echo ""
echo "------------------------------------------------------------"
