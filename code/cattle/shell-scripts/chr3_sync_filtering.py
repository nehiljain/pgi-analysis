
import sys
import pandas as pd
import numpy as np
import csv
from os import listdir
from os.path import isfile, join

chr_no = '3'
snplist_filepath = r'/share/volatile_scratch/nehil/mouse/filtered/snp_list_chr'+ chr_no +'.csv'
print 'reading the snp list in memory'
chr_snp_list = pd.read_csv(snplist_filepath)
print chr_snp_list.shape
sync_dir_path = r'/share/volatile_scratch/nehil/pgi_wc/mouse/dna_seq/genomequebec_cmh_sync/'
sync_filepath = sync_dir_path + 'allSamples.chr3:1-160039680.sync'

out_file_name = r'/share/volatile_scratch/nehil/mouse/filtered_sync/filtered_' + 'chr_' + chr_no + '.sync'
print out_file_name + '\n'

if chr_snp_list.shape[0] != 0:
  pos_list = chr_snp_list.pos.values
  snp_max_pos = np.amax(pos_list)
  print snp_max_pos
  with open(sync_filepath,'rb') as tsvin, open(out_file_name, 'wb') as csvout:
    tsvin = csv.reader(tsvin, delimiter='\t')
    csvout = csv.writer(csvout, delimiter='\t')
    for row in tsvin:
        snp_pos = int(row[1])
        ar = np.where(pos_list == snp_pos)[0]
        if snp_pos > snp_max_pos:
            print "breaking the read loop" + snp_pos
            break
        print snp_pos
        if ar.size != 0:
            if ar[0].size != 0:
              csvout.writerow(row)
              print row
                           13,1          All
