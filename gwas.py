import os

fileObject = open("cmh-to-gwas.sh", "wb")



fileObject.write("#!/bin/sh \n")

for file in os.listdir("cmh-analysis-files"):
        if file.endswith(".cmh"):
                fileObject.write("perl ~/popoolation2_1201/export/cmh2gwas.pl"+
                 " --input cmh-analysis-files/" + file + " --output cmh-igv-files/" + os.path.splitext(file)[0] + ".gwas --min-pvalue 0 \n\n")
                	
fileObject.close()