#__author__ = 'nehiljain'


"""Functions of all the sheel commands in the pipeline which will be used as tasks in Ruffus.
"""

import sys
import subprocess
import os
import ruffus
import re

from ruffus import *

PICARD_JAR = "/home/kasia/picard/picard.jar"
POPOOLATION2_DIR = "/home/kzukowski/soft/popoolation2_1201/"
CATTLE_REF_DIR = "/home/kasia/data/reference/"
CATTLE_REF_SEQ_FILE = CATTLE_REF_DIR + "Bos_taurus.UMD3.1.dna.toplevel.fa"

INIT_DIR = "/home/kasia/data/gq_alignment_bam/pooled-data"
LOG_DIR = "/home/kasia/data/logs/"

# BASE_DIR = "/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/"

# CLEANSAM_OUT_DIR =  "nehil_cleansam_bam"
# MAPQ20_OUT_DIR =  "nehil_mapq20_bam"
# FIXMATE_OUT_DIR =  "nehil_fixmate_bam"
# DEDUP_OUT_DIR =  "nehil_mark_duplicates_bam"
# MULTIPLE_METRICS_OUT_DIR =  """nehil_collect_multiple_metrics_CollectMultipleMetrics/"""
# COLLECT_GC_BIAS_METRICS_OUT_DIR =  """ nehil_collect_gc_bias_meterics_CollectGcBiasMetrics/"""


# SYNC_OUT_DIR =  "nehil_samtools_mpileup/"
# SYNC_OUT_DIR =  "nehil_mpileup2sync_sync/"


def ensure_path_exists(path):
    """If the directory does not exist it creates one.
    """
    if not os.path.exists(path):
        os.makedirs(path)


def run_cmd(cmd_str, output_log_file = 'output.log',
            error_log_file = 'error.log'):
    """Runs the command as given in command string.

    This function uses subprocess to run shell commands in cmd_str. Throws an exception if run
    command fails.
    returns stdout and stderror in a list"""
    print "-"*100
    print "Running Shell Command: \n"
    print cmd_str
    o_log = open(output_log_file, 'w')
    e_log = open(error_log_file, 'w')
    process_returncode = subprocess.call(cmd_str, stdout = o_log,
                               stderr = e_log, shell = True)
    print "FINISHED"
    print process_returncode
    print "-"*100

    o_log.close()
    e_log.close()
    if process_returncode != 0:
        raise Exception("Failed to run '%s' \n Non-zero exit status %s" %
                        (cmd_str, process_returncode))


def get_all_init_filepaths(dir_path):
    """This function gets the BAM files from the INIT DIR folder.
    """
    pattern = re.compile("Fixmate.bam");
    init_bam_files = [ f for f in os.listdir(dir_path)
                        if re.search(pattern, f)]
    print "-"*100
    print("INIT Files")
    print(init_bam_files)
    print "-"*100
    return init_bam_files

init_files = get_all_init_filepaths(INIT_DIR)
print("init_files",init_files)



@transform(init_files,
           suffix(".bam"),
           ".CleanSam.bam")
def picard_cleansam(init_files, output_file_name):
    """
    :param input_file: string '12766.sorted.bam'
    :param output_file: string '12766.CleanSam.bam'
    :return: stderror and stdout
    """
    print("INIT Files", init_files)
    out_log_file_path = LOG_DIR + os.path.splitext(init_files)[0] + ".out.log"
    err_log_file_path = LOG_DIR + os.path.splitext(init_files)[0] + ".err.log"
    out_file_path = output_file_name

    command_str = ("""java -Xmx8g -jar {picard} CleanSam INPUT={inp} OUTPUT={outp} VALIDATION_STRINGENCY=SILENT """
<<<<<<< HEAD
    """ CREATE_INDEX=true TMP_DIR=/tmp""".format(picard = PICARD_JAR, inp = init_files, outp = out_file_path))
=======
        """ CREATE_INDEX=true TMP_DIR=/tmp""".format(picard = PICARD_JAR,
        inp = init_files<script src="/bower_components/readmore/readmore.min.js"></script>, outp = out_file_path))
>>>>>>> 7e4dd3e577757d9505b032099369d65402240787
    run_cmd(command_str, out_log_file_path, err_log_file_path)


# @follows(picard_cleansam)
# @transform(picard_cleansam,
#            suffix(".CleanSam.bam"),
#            ".MAPQ20.bam")
# def samtools_mapq20(input_file, output_file):
#     """
#     :param input_file:
#     :param output_file:
#     :param file_name:
#     :return:
#     """
#     out_log_file_path = LOG_DIR + os.path.splitext(input_file)[0] + ".out.log"
#     err_log_file_path = LOG_DIR + os.path.splitext(input_file)[0] + ".err.log"

#     command_str = ("""samtools view -bq 20 {inp}""".format(
#       inp = input_file))
#     run_cmd(command_str, output_file, err_log_file_path)


# @follows(samtools_mapq20)
# @transform(init_files,
#            suffix(".MAPQ20.bam"),
#            ".picard.Fixmate.bam")
# def picard_fixmate(input_file, output_file):
#     """

#     :param input_file:
#     :param output_file:
#     :param log_file:
#     :return:
#     """
#     out_log_file_path = LOG_DIR + os.path.splitext(input_file)[0] + ".out.log"
#     err_log_file_path = LOG_DIR + os.path.splitext(input_file)[0] + ".err.log"

#     command_str = ("""java -Xmx8g -jar {picard} FixMateInformation INPUT={inp} OUTPUT={outp} VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate ASSUME_SORTED=true ADD_MATE_CIGAR=true """
#         """CREATE_INDEX=true TMP_DIR=/tmp""".format(picard = PICARD_JAR,
#         inp = input_file, outp = output_file))

#     run_cmd(command_str, out_log_file_path, err_log_file_path)

# @follows(picard_fixmate)
@transform(init_files,
           suffix(".picard.Fixmate.bam"),
           ".picard.DeDup.bam")

def picard_mark_duplicates(init_file, output_file):
    """

    :param input_file:
    :param output_file:
    :return:
    """
    print("INIT Files", init_files, init_file)
    out_metrics_file_path = os.path.splitext(init_file)[0] + ".metrics"
    out_log_file_path = LOG_DIR + os.path.splitext(init_file)[0] + ".out.log"
    err_log_file_path = LOG_DIR + os.path.splitext(init_file)[0] + ".err.log"

    command_str = ("""java -jar {picard} MarkDuplicates INPUT={inp} OUTPUT={outp} METRICS_FILE={metrp} VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true ASSUME_SORTED=true """
        """ CREATE_INDEX=true TMP_DIR=/tmp""".format(picard = PICARD_JAR,
        inp = init_file, outp = output_file,
        metrp = out_metrics_file_path ))
    run_cmd(command_str, out_log_file_path, err_log_file_path)

# def picard_collect_multiple_metrics(input_file, output_file, file_name):
#     """

#     :param input_file:
#     :param output_file:
#     :param log_file:
#     :return:
#     """
#     out_file_path = MULTIPLE_METRICS_OUT_DIR + "12429.picard" + ".CollectMultipleMetrics"
#     out_log_file_path = LOG_DIR + "12429.picard.CollectMultipleMetrics" + ".out.log"
#     err_log_file_path = LOG_DIR + "12429.picard.CollectMultipleMetrics" + ".err.log"
#     in_file_path = DEDUP_OUT_DIR + "12429.picard.DeDup" + ".bam"

#     command_str = ("""java -Xmx8g -jar {picard} CollectMultipleMetrics INPUT={inp} OUTPUT={outp} VALIDATION_STRINGENCY=SILENT TMP_DIR=/tmp""".format(picard = PICARD_JAR, inp = in_file_path, outp = out_file_path))
#     print command_str


# def picard_collect_gc_bias_metrics(input_file, output_file, file_name):
#     """

#     :param input_file:
#     :param output_file:
#     :param log_file:
#     :return:
#     """
#     out_file_path = (COLLECT_GC_BIAS_METRICS_OUT_DIR
#         + "12429.picard" +    ".CollectGcBiasMetrics")
#     out_chart_file_path = (COLLECT_GC_BIAS_METRICS_OUT_DIR
#     + "12429.picard" + ".CollectGcBiasMetrics" + ".pdf")
#     out_log_file_path = (LOG_DIR + "12429.picard.CollectGcBiasMetrics"
#     +        ".out.log")
#     err_log_file_path = (LOG_DIR + "12429.picard.CollectGcBiasMetrics"
#     +   ".err.log")
#     in_file_path = DEDUP_OUT_DIR + "12429.picard.DeDup" + ".bam"

#     command_str = ("""java -Xmx8g -jar {picard} CollectMultipleMetrics INPUT={inp} OUTPUT={outp} CHART_OUTPUT={chartp} """
#     """ REFERENCE_SEQUENCE={ref_seq} """
#     """ VALIDATION_STRINGENCY=SILENT TMP_DIR=/tmp""".format(picard = PICARD_JAR, inp = in_file_path, outp = out_file_path, ref_seq = CATTLE_REF_SEQ_FILE ,
#          chartp = out_chart_file_path))
#     print command_str



# def samtools_mpileup(input_files,
#                               output_file,
#                               reference_seq,
#                               log_file):
#     """

#     :param input_file:
#     :param output_file:
#     :param log_file:
#     :return:
#     """
#     pass


# def popoolation2_mpileup_to_sync(input_file, output_file, log_file):
#     """
#     :param input_file:
#     :param output_file:
#     :param log_file:
#     :return:
#     """
#     mpileup_path = POPOOLATION2_DIR + "mpileup2sync.jar"
#     in_file_path = ("/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.mpileup")
#     log_file_path = "/share/volatile_scratch/kzukowski/pgi/cattle/test/log/12766.mpileup2sync.log"
#     out_file_path = "/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.sync"

#     command_str = ("""java -Xmx8g -jar {mpileup} --input {inp}
# --output {outp} --fastq-type sanger --min-qual 20 --threads 2"""
#         """CREATE_INDEX=true TMP_DIR=/tmp""".format(mpileup = mpileup_path, inp = in_file_path, outp = out_file_path))

#     logs = run_cmd(command_str)
#     #write to log file
#     print(logs)
#     pass




ensure_path_exists(LOG_DIR)

os.chdir(INIT_DIR)
pipeline_get_task_names()
# pipeline_printout(output_stream = sys.stdout, target_tasks = [ picard_mark_duplicates], verbose=10,checksum_level=3, verbose_abbreviated_path=1)
pipeline_run(target_tasks = [ picard_mark_duplicates], verbose=6,checksum_level=3,verbose_abbreviated_path=1)






