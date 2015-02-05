#!/bin/bash
clear

echo ""
echo "Host: `hostname`"
echo "Login: $USER"
echo "Run time: `date`"
echo ""

date
echo
echo "popoolation2/cmh-test.pl"
echo

perl /home/kzukowski/soft/popoolation2_1201/cmh-test.pl \
--input /share/volatile_scratch/nehil/pgi_wc/mouse/dna_seq/nehil_filtered_sync/filtered_chr_3.sync \
--output /share/volatile_scratch/kzukowski/pgi/mouse/test/chr3.cmh \
--min-count -1 \
--min-coverage -1 \
--max-coverage 1000000000 \
--remove-temp \
--population 3-1,7-2,6-8,5-4 \
>/share/volatile_scratch/kzukowski/pgi/mouse/test/log/12766.mpileup2sync.log 2>&1

cut -f 1,2,3,12 /share/volatile_scratch/kzukowski/pgi/mouse/test/chr3.cmh >/share/volatile_scratch/kzukowski/pgi/mouse/test/chr3.gwas

################ REMOVING FILES -r !!!
################ NEHIL STEADY WITH THIS function!! :D
### not used due genotype distribution plots?

#rm -r /share/volatile_scratch/kzukowski/pgi/mouse/test/chr3.cmh
#rm -r /share/volatile_scratch/kzukowski/pgi/mouse/test/chr3.cmh.params

