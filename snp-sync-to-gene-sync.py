import os

fileObject = open("~/pgi-ngs-analysis/snp-to-gene-sync.sh", "wb")



fileObject.write("#!/bin/sh \n")

for file in os.listdir("/share/volatile_scratch/nehil/sync_files"):
        if file.endswith(".sync"):
        		# fileObject.write("echo \" "+ str(file) + "\"\n\n") 
                fileObject.write("perl -d:DProf ~/popoolation2_1201/create-genewise-sync.pl " + 
                "--input /share/volatile_scratch/nehil/sync_files/" + file + 
                " --gtf /share/volatile_scratch/nehil/sync_files/Mus_musculus.GRCm38.75.gtf" + 
                " --output /share/volatile_scratch/nehil/gene-sync-files/" + os.path.splitext(file)[0] + ".sync  \n\n")

fileObject.close()