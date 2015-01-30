#__author__ = 'nehiljain'


"""Functions of all the sheel commands in the pipeline which will be used as tasks in Ruffus.
"""

import sys
import subprocess
import os
import ruffus

from ruffus import *

PICARD_JAR = "/share/apps/picard/git/bin/picard.jar"
POPOOLATION2_DIR = "/home/kzukowski/soft/popoolation2_1201/"
CATTLE_REF_DIR = "/share/volatile_scratch/kzukowski/pgi/cattle/reference/"
CATTLE_REF_SEQ_FILE = CATTLE_REF_DIR + "Bos_taurus.UMD3.1.dna.toplevel.fa"

INIT_DIR = "/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/gq_alignment_bam/"

BASE_DIR = "/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/"
LOG_DIR = "/share/volatile_scratch/nehil/pgi_wc/logs/"
CLEANSAM_OUT_DIR = BASE_DIR + "nehil_cleansam_bam/"
MAPQ20_OUT_DIR = BASE_DIR + "nehil_mapq20_bam/"
FIXMATE_OUT_DIR = BASE_DIR + "nehil_fixmate_bam/"
DEDUP_OUT_DIR = BASE_DIR + "nehil_mark_duplicates_bam/"
MULTIPLE_METRICS_OUT_DIR = BASE_DIR + """nehil_collect_multiple_metrics_CollectMultipleMetrics/"""
COLLECT_GC_BIAS_METRICS_OUT_DIR = BASE_DIR + """ nehil_collect_gc_bias_meterics_CollectGcBiasMetrics/"""


SYNC_OUT_DIR = BASE_DIR + "nehil_samtools_mpileup/"
SYNC_OUT_DIR = BASE_DIR + "nehil_mpileup2sync_sync/"


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
    init_bam_files = [ f for f in os.listdir(dir_path)
                        if os.path.splitext(f)[1] == '.bam']
    return init_bam_files

init_files = []

def init_stub():
    pass


@follows(init_stub, mkdir(CLEANSAM_OUT_DIR))
@transform(init_files, suffix(".bam"),
           [".CleanSam.bam", ".CleanSam.out.log", ".CleanSam.err.log"])
def picard_cleansam(input_file, output_file_names):
    """
    :param input_file: string '12766.sorted.bam'
    :param output_file: string '12766.CleanSam.bam'
    :return: stderror and stdout
    """

    in_file_path = INIT_DIR + init_files
    out_log_file_path = LOG_DIR + output_file_names[1]
    err_log_file_path = LOG_DIR + output_file_names[2]
    out_file_path = CLEANSAM_OUT_DIR + output_file_names[0]

    command_str = ("""java -Xmx8g -jar {picard} CleanSam INPUT={inp} OUTPUT={outp} VALIDATION_STRINGENCY=SILENT """
        """ CREATE_INDEX=true TMP_DIR=/tmp""".format(picard = PICARD_JAR,
        inp = in_file_path, outp = out_file_path))
    run_cmd(command_str, out_log_file_path, err_log_file_path)


@follows("picard_cleansam", mkdir(MAPQ20_OUT_DIR))
@transform(picard_cleansam, suffix(".bam"),
           [".MAPQ20.bam", ".MAPQ20.out.log", ".MAPQ20.err.log"])
def samtools_mapq20(input_file, output_file, filename):
    """

    :param input_file:
    :param output_file:
    :param file_name:
    :return:
    """
    out_log_file_path = LOG_DIR + output_file_names[1]
    err_log_file_path = LOG_DIR + output_file_names[2]
    out_file_path = MAPQ20_OUT_DIR + output_file_names[0]

    command_str = ("""samtools view -bq 20 {inp}""".format(
        inp = input_file))
    run_cmd(command_str, out_file_path, err_log_file_path)


def picard_fixmate(input_file, output_file, filename):
    """

    :param input_file:
    :param output_file:
    :param log_file:
    :return:
    """
    in_file_path = MAPQ20_OUT_DIR + "12429.samtools.MAPQ20" + ".bam"
    out_log_file_path = LOG_DIR + "12429.picard.Fixmate" + ".out.log"
    err_log_file_path = LOG_DIR + "12429.picard.Fixmate" + ".err.log"
    out_file_path = FIXMATE_OUT_DIR + "12429.picard.Fixmate" + ".bam"

    command_str = ("""java -Xmx8g -jar {picard} FixMateInformation INPUT={inp} OUTPUT={outp} VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate ASSUME_SORTED=true ADD_MATE_CIGAR=true """
        """CREATE_INDEX=true TMP_DIR=/tmp""".format(picard = PICARD_JAR,
        inp = in_file_path, outp = out_file_path))
    print command_str
    # run_cmd(command_str, out_log_file_path, err_log_file_path)



def picard_mark_duplicates(input_file, output_file, filename):
    """

    :param input_file:
    :param output_file:
    :param log_file:
    :return:
    """
    out_file_path = DEDUP_OUT_DIR + "12429.picard.DeDup" + ".bam"
    out_metrics_file_path = DEDUP_OUT_DIR + "12429.picard.DeDup" + ".metrics"
    out_log_file_path = LOG_DIR + "12429.picard.DeDup" + ".out.log"
    err_log_file_path = LOG_DIR + "12429.picard.DeDup" + ".err.log"
    in_file_path = FIXMATE_OUT_DIR + "12429.picard.Fixmate" + ".bam"

    command_str = ("""java -Xmx8g -jar {picard} MarkDuplicates INPUT={inp} OUTPUT={outp} METRICS_FILE={metrp} VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true ASSUME_SORTED=true """
        """ CREATE_INDEX=true TMP_DIR=/tmp""".format(picard = PICARD_JAR,
        inp = in_file_path, outp = out_file_path,
        metrp = out_metrics_file_path ))
    print command_str


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





if __name__ == '__main__':

    init_files = get_all_init_filepaths(INIT_DIR)
    ensure_path_exists(BASE_DIR)
    ensure_path_exists(LOG_DIR)
    ensure_path_exists(CLEANSAM_OUT_DIR)
    ensure_path_exists(MAPQ20_OUT_DIR)
    ensure_path_exists(DEDUP_OUT_DIR)
    ensure_path_exists(FIXMATE_OUT_DIR)
    ensure_path_exists(MULTIPLE_METRICS_OUT_DIR)
    ensure_path_exists(COLLECT_GC_BIAS_METRICS_OUT_DIR)
    init_files = get_all_init_filepaths(INIT_DIR)
    pipeline_run(target_tasks = [samtools_mapq20])
    # picard_cleansam(1,2,3)
    # samtools_mapq20(1,2,3)
    # picard_fixmate(1,2,3)
    # picard_mark_duplicates(1,2,3)
    # picard_collect_multiple_metrics(1,2,3)
    # picard_collect_gc_bias_metrics(1,2,3)
    sys.exit(0)






