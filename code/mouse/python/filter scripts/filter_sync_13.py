
import sys
import pandas as pd
import numpy as np
import csv
from os import listdir
from os.path import isfile, join

chr_no = '13'
snplist_filepath = r'/share/volatile_scratch/nehil/pgi_wc/mouse/dna_seq/kacper_filtered_snp_list_csv'+ chr_no +'.csv'

ensure_path_exists(snplist_filepath)
print 'reading the snp list in memory'
chr_snp_list = pd.read_csv(snplist_filepath)
print chr_snp_list.shape
sync_dir_path = r'/share/volatile_scratch/nehil/pgi_wc/mouse/dna_seq/gq_cmh_sync'
ensure_path_exists(sync_dir_path)
sync_filepath = sync_dir_path + 'allSamples.chr13:1-120421639.sync'

out_file_name = r'/share/volatile_scratch/nehil/pgi_wc/mouse/dna_seq/nehil_filtered_sync' + 'chr_' + chr_no + '.sync'
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
            print "breaking the read loop" + str(snp_pos)
            break
        print snp_pos
        if ar.size != 0:
            if ar[0].size != 0:
              csvout.writerow(row)
              print row



def ensure_path_exists(path):
    """If the directory does not exist it creates one.
    """
    if not os.path.exists(path):
        os.makedirs(path)
