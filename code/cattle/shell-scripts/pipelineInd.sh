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
beagle="java -jar /home/kzukowski/soft/beagle/beagle.r1399.jar"
echo "beagle"
beagle3="java -jar /home/kzukowski/soft/beagle/beagle3.jar"
echo "beagle3"
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
echo ""

echo "samples list:"
cat $path_vcf/sample.list
echo ""

echo "============================================================"
echo "analysis:"

for s in `cat $path_vcf/sample.list`;
	do
	echo $s >$path_vcf/sample

	begin=$(date +%s)
	$picard SortVcf \
	I=$path_vcf0/$s.hc.vcf \
	O=$path_vcf/$s.hc.sort.vcf \
	SEQUENCE_DICTIONARY=/share/volatile_scratch/nehil/pgi_wc/cattle/reference/77/Bos_taurus.UMD3.1.dna.toplevel.dict \
	>$path_log/$s.picard.SortVcf.log 2>&1
	termin=$(date +%s)
	difftime=$(( $termin - $begin ))
	echo "picard SortVcf sample $s runtime: $difftime seconds"

	done
	rm -r $path_vcf/sample

echo ""
echo "------------------------------------------------------------"

begin=$(date +%s)
$gatk -T CombineVariants -R $path_ref \
-V:12429 $path_vcf/12429.hc.sort.vcf \
-V:12569 $path_vcf/12569.hc.sort.vcf \
-V:12745 $path_vcf/12745.hc.sort.vcf \
-V:12766 $path_vcf/12766.hc.sort.vcf \
-V:24237 $path_vcf/24237.hc.sort.vcf \
-V:24406 $path_vcf/24406.hc.sort.vcf \
-V:32714 $path_vcf/32714.hc.sort.vcf \
-V:51844 $path_vcf/51844.hc.sort.vcf \
-V:51912 $path_vcf/51912.hc.sort.vcf \
-V:53763 $path_vcf/53763.hc.sort.vcf \
-V:59018 $path_vcf/59018.hc.sort.vcf \
-V:60070 $path_vcf/60070.hc.sort.vcf \
-V:60113 $path_vcf/60113.hc.sort.vcf \
-V:60383 $path_vcf/60383.hc.sort.vcf \
-V:61335 $path_vcf/61335.hc.sort.vcf \
-V:61536 $path_vcf/61536.hc.sort.vcf \
-V:62206 $path_vcf/62206.hc.sort.vcf \
-V:62225 $path_vcf/62225.hc.sort.vcf \
-V:62227 $path_vcf/62227.hc.sort.vcf \
-V:62242 $path_vcf/62242.hc.sort.vcf \
-V:62668 $path_vcf/62668.hc.sort.vcf \
-V:62696 $path_vcf/62696.hc.sort.vcf \
-V:62727 $path_vcf/62727.hc.sort.vcf \
-V:62770 $path_vcf/62770.hc.sort.vcf \
-V:62797 $path_vcf/62797.hc.sort.vcf \
-V:63382 $path_vcf/63382.hc.sort.vcf \
-V:157222 $path_vcf/157222.hc.sort.vcf \
-V:157240 $path_vcf/157240.hc.sort.vcf \
-V:157243 $path_vcf/157243.hc.sort.vcf \
-V:157250 $path_vcf/157250.hc.sort.vcf \
-V:157254 $path_vcf/157254.hc.sort.vcf \
-V:157257 $path_vcf/157257.hc.sort.vcf \
-V:157258 $path_vcf/157258.hc.sort.vcf \
-V:157269 $path_vcf/157269.hc.sort.vcf \
-V:157275 $path_vcf/157275.hc.sort.vcf \
-V:157278 $path_vcf/157278.hc.sort.vcf \
-V:157287 $path_vcf/157287.hc.sort.vcf \
-V:157289 $path_vcf/157289.hc.sort.vcf \
-V:157291 $path_vcf/157291.hc.sort.vcf \
-V:157292 $path_vcf/157292.hc.sort.vcf \
-V:157295 $path_vcf/157295.hc.sort.vcf \
-V:157305 $path_vcf/157305.hc.sort.vcf \
-V:157309 $path_vcf/157309.hc.sort.vcf \
-V:157316 $path_vcf/157316.hc.sort.vcf \
-V:157337 $path_vcf/157337.hc.sort.vcf \
-V:157341 $path_vcf/157341.hc.sort.vcf \
-V:157345 $path_vcf/157345.hc.sort.vcf \
-V:157349 $path_vcf/157349.hc.sort.vcf \
-V:157354 $path_vcf/157354.hc.sort.vcf \
-V:157357 $path_vcf/157357.hc.sort.vcf \
-V:157362 $path_vcf/157362.hc.sort.vcf \
-V:157368 $path_vcf/157368.hc.sort.vcf \
-V:157374 $path_vcf/157374.hc.sort.vcf \
-V:157375 $path_vcf/157375.hc.sort.vcf \
-V:157376 $path_vcf/157376.hc.sort.vcf \
-V:157394 $path_vcf/157394.hc.sort.vcf \
-V:157399 $path_vcf/157399.hc.sort.vcf \
-V:157402 $path_vcf/157402.hc.sort.vcf \
-V:157981 $path_vcf/157981.hc.sort.vcf \
-V:157985 $path_vcf/157985.hc.sort.vcf \
-V:157987 $path_vcf/157987.hc.sort.vcf \
-V:12441 $path_vcf/12441.hc.sort.vcf \
-V:12718 $path_vcf/12718.hc.sort.vcf \
-V:23934 $path_vcf/23934.hc.sort.vcf \
-V:24321 $path_vcf/24321.hc.sort.vcf \
-V:24403 $path_vcf/24403.hc.sort.vcf \
-V:32717 $path_vcf/32717.hc.sort.vcf \
-V:32721 $path_vcf/32721.hc.sort.vcf \
-V:51822 $path_vcf/51822.hc.sort.vcf \
-V:51866 $path_vcf/51866.hc.sort.vcf \
-V:51874 $path_vcf/51874.hc.sort.vcf \
-V:53764 $path_vcf/53764.hc.sort.vcf \
-V:58635 $path_vcf/58635.hc.sort.vcf \
-V:60130 $path_vcf/60130.hc.sort.vcf \
-V:61334 $path_vcf/61334.hc.sort.vcf \
-V:62103 $path_vcf/62103.hc.sort.vcf \
-V:62139 $path_vcf/62139.hc.sort.vcf \
-V:62156 $path_vcf/62156.hc.sort.vcf \
-V:62196 $path_vcf/62196.hc.sort.vcf \
-V:62218 $path_vcf/62218.hc.sort.vcf \
-V:62243 $path_vcf/62243.hc.sort.vcf \
-V:62254 $path_vcf/62254.hc.sort.vcf \
-V:62312 $path_vcf/62312.hc.sort.vcf \
-V:62666 $path_vcf/62666.hc.sort.vcf \
-V:62703 $path_vcf/62703.hc.sort.vcf \
-V:62761 $path_vcf/62761.hc.sort.vcf \
-V:62866 $path_vcf/62866.hc.sort.vcf \
-V:62868 $path_vcf/62868.hc.sort.vcf \
-V:63193 $path_vcf/63193.hc.sort.vcf \
-V:63247 $path_vcf/63247.hc.sort.vcf \
-V:63624 $path_vcf/63624.hc.sort.vcf \
-V:63641 $path_vcf/63641.hc.sort.vcf \
-V:69299 $path_vcf/69299.hc.sort.vcf \
-V:157214 $path_vcf/157214.hc.sort.vcf \
-V:157221 $path_vcf/157221.hc.sort.vcf \
-V:157223 $path_vcf/157223.hc.sort.vcf \
-V:157228 $path_vcf/157228.hc.sort.vcf \
-V:157229 $path_vcf/157229.hc.sort.vcf \
-V:157282 $path_vcf/157282.hc.sort.vcf \
-V:157303 $path_vcf/157303.hc.sort.vcf \
-V:157304 $path_vcf/157304.hc.sort.vcf \
-V:157306 $path_vcf/157306.hc.sort.vcf \
-V:157310 $path_vcf/157310.hc.sort.vcf \
-V:157318 $path_vcf/157318.hc.sort.vcf \
-V:157320 $path_vcf/157320.hc.sort.vcf \
-V:157324 $path_vcf/157324.hc.sort.vcf \
-V:157328 $path_vcf/157328.hc.sort.vcf \
-V:157329 $path_vcf/157329.hc.sort.vcf \
-V:157342 $path_vcf/157342.hc.sort.vcf \
-V:157348 $path_vcf/157348.hc.sort.vcf \
-V:157355 $path_vcf/157355.hc.sort.vcf \
-V:157379 $path_vcf/157379.hc.sort.vcf \
-V:157395 $path_vcf/157395.hc.sort.vcf \
-V:157397 $path_vcf/157397.hc.sort.vcf \
-V:157989 $path_vcf/157989.hc.sort.vcf \
-V:157994 $path_vcf/157994.hc.sort.vcf \
-V:157996 $path_vcf/157996.hc.sort.vcf \
-V:158002 $path_vcf/158002.hc.sort.vcf \
-V:158003 $path_vcf/158003.hc.sort.vcf \
-V:158009 $path_vcf/158009.hc.sort.vcf \
-V:158015 $path_vcf/158015.hc.sort.vcf \
-V:158017 $path_vcf/158017.hc.sort.vcf \
-V:158018 $path_vcf/158018.hc.sort.vcf \
-V:158028 $path_vcf/158028.hc.sort.vcf \
-V:158036 $path_vcf/158036.hc.sort.vcf \
-V:158038 $path_vcf/158038.hc.sort.vcf \
-genotypeMergeOptions UNIQUIFY \
-o $path_vcfC/combined.vcf \
>$path_log/gatk.CombineVariants.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk CombineVariants runtime: $difftime seconds"

begin=$(date +%s)
$gatk -T VariantRecalibrator -R $path_ref \
-input $path_vcfC/combined.vcf \
-resource:VCF,known=true,training=true,truth=true,prior=2.0 /share/volatile_scratch/nehil/pgi_wc/cattle/reference/77/snp.Bos_taurus.vcf \
-an DP \
-an QD \
-an FS \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile $path_stats/recalibrate_SNP.recal \
-tranchesFile $path_stats/recalibrate_SNP.tranches \
-rscriptFile $path_stats/recalibrate_SNP_plots.R \
>$path_log/gatk.VariantRecalibratorSNP.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk VariantRecalibrator SNP runtime: $difftime seconds"

begin=$(date +%s)
$gatk -T ApplyRecalibration -R $path_ref \
-input $path_vcfC/combined.vcf \
-mode SNP \
--ts_filter_level 99.0 \
-recalFile $path_stats/recalibrate_SNP.recal \
-tranchesFile $path_stats/recalibrate_SNP.tranches \
-o $path_vcfC/recalibrated_snps_raw_indels.vcf \
>$path_log/gatk.ApplyRecalibrationSNP.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk ApplyRecalibration SNP runtime: $difftime seconds"

begin=$(date +%s)
$gatk -T VariantRecalibrator -R $path_ref \
-input $path_vcfC/recalibrated_snps_raw_indels.vcf \
-resource:VCF,known=true,training=true,truth=true,prior=12.0 /share/volatile_scratch/nehil/pgi_wc/cattle/reference/77/indels.Bos_taurus.vcf \
-an QD \
-an DP \
-an FS \
-an MQRankSum \
-an ReadPosRankSum \
-mode INDEL \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
--maxGaussians 4 \
-recalFile $path_stats/recalibrate_INDEL.recal \
-tranchesFile $path_stats/recalibrate_INDEL.tranches \
-rscriptFile $path_stats/recalibrate_INDEL_plots.R \
>$path_log/gatk.VariantRecalibratorINDEL.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk VariantRecalibrator INDEL runtime: $difftime seconds"

begin=$(date +%s)
$gatk -T ApplyRecalibration -R $path_ref \
-input $path_vcfC/recalibrated_snps_raw_indels.vcf \
-mode INDEL \
--ts_filter_level 99.0 \
-recalFile $path_stats/recalibrate_INDEL.recal \
-tranchesFile $path_stats/recalibrate_INDEL.tranches \
-o $path_vcfC/recalibrated_variants.vcf \
>$path_log/gatk.ApplyRecalibrationINDEL.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk ApplyRecalibration INDEL runtime: $difftime seconds"

begin=$(date +%s)
$gatk -T SelectVariants -R $path_ref \
-V $path_vcfC/recalibrated_variants.vcf \
-selectType SNP \
-o $path_vcfC/raw_snps.vcf \
>$path_log/gatk.SelectVariantsSNP.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk SelectVariants SNP runtime: $difftime seconds" 

begin=$(date +%s)
$gatk -T VariantFiltration -R $path_ref \
-V $path_vcfC/raw_snps.vcf \
-filterName FS -filter "FS > 60.0" \
-filterName QD -filter "QD < 2.0" \
-filterName MQ -filter "MQ < 40.0" \
-filterName HaplotypeScore -filter "HaplotypeScore > 13.0" \
-filterName MappingQualityRankSum -filter "MappingQualityRankSum < -12.5" \
-filterName ReadPosRankSum -filter "ReadPosRankSum -8.0" \
-filterName QUAL -filter "QUAL < 30" \
-o $path_vcfC/filtered_snps.vcf \
>$path_log/gatk.VariantFiltrationSNP.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk VariantFiltration SNP runtime: $difftime seconds"

begin=$(date +%s)
$gatk -T SelectVariants -R $path_ref \
-V $path_vcfC/recalibrated_variants.vcf \
-selectType INDEL \
-o $path_vcfC/raw_indels.vcf \
>$path_log/gatk.SelectVariantsINDEL.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk SelectVariants INDEL runtime: $difftime seconds"

begin=$(date +%s)
$gatk -T VariantFiltration -R $path_ref \
-V $path_vcfC/raw_indels.vcf \
-filterName QD -filter "QD < 2.0" \
-filterName FS -filter "FS > 200.0" \
-filterName ReadPosRankSum -filter "ReadPosRankSum < -20.0" \
-o $path_vcfC/filtered_indels.vcf \
>$path_log/gatk.VariantFiltrationINDEL.log 2>&1
termin=$(date +%s)
difftime=$(( $termin - $begin ))
echo "gatk VariantFiltration INDEL runtime: $difftime seconds"

