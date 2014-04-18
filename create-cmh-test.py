import os

fileObject = open("cmh-test-all.sh", "wb")



fileObject.write("#!/bin/sh \n")

for file in os.listdir("/share/volatile_scratch/nehil/sync_files"):
        if file.endswith(".sync"):
                fileObject.write("perl -d:DProf ~/popoolation2_1201/cmh-test.pl " + 
                	"--input /share/volatile_scratch/nehil/sync_files/" + file +" --output cmh-analysis-files/" + os.path.splitext(file)[0] + ".cmh " + 
                	"--min-count 12 --min-coverage 20 --max-coverage 10000000 --population 6-1,6-2,6-4,6-8,3-1,3-2,3-8,3-4,7-1,7-8,7-2,7-4,5-1,5-2,5-8,5-4 \n\n" )
fileObject.close()

