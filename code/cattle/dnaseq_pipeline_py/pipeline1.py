#__author__ = 'nehiljain'


"""Functions of all the sheel commands in the pipeline which will be used as tasks in Ruffus.
"""

import sys
import subprocess
import os
#import ruffus

PICARD_JAR = "/share/apps/picard/git/bin/picard.jar"
POPOOLATION2_DIR = "/home/kzukowski/soft/popoolation2_1201/"
CATTLE_REF_DIR = "/share/volatile_scratch/kzukowski/pgi/cattle/reference/"
CATTLE_REF_SEQ_FILE = CATTLE_REF_DIR + "Bos_taurus.UMD3.1.dna.toplevel.fa"

INIT_DIR = "/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/gq_alignment_bam/"

BASE_DIR = "/share/volatile_scratch/nehil/pgi_wc/cattle/dna_seq/"
LOG_DIR = "/share/volatile_scratch/nehil/pgi_wc/logs/"
CLEANSAM_OUT_DIR = BASE_DIR + "nehil_cleansam_bam/"
MAPQ20_OUT_DIR =   BASE_DIR + "nehil_mapq20_bam/"

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
    print "Running Shell Command: \n\n"
    print "cmd_str"
    o_log = open(output_log_file, 'w')
    e_log = open(error_log_file, 'w')
    process_returncode = subprocess.call(cmd_str, stdout = output_log_file,
                               stderr = error_log_file, shell = True)
    print "FINISHED"
    print process_returncode
    print "-"*100

    o_log.close()
    e_log.close()
    if process_returncode != 0:
        raise Exception("Failed to run '%s' \n Non-zero exit status %s" %
                        (cmd_str, process_returncode))


def picard_cleansam(input_file, output_file, file_name):
    """
    :param input_file: string '12766.sorted.bam'
    :param output_file: string '12766.CleanSam.bam'
    :param log_file: string '12766.picard.CleanSam.log'
    :return: stderror and stdout
    """

    in_file_path = BASE_DIR + "gq_alignment_bam/12429.sorted.bam"
    out_log_file_path = LOG_DIR + "12429.picard.CleanSam" + ".out.log"
    err_log_file_path = LOG_DIR + "12429.picard.CleanSam" + ".err.log"
    out_file_path = CLEANSAM_OUT_DIR + "/12429.picard.CleanSam.bam"

    command_str = ("""java -jar {picard} CleanSam INPUT={inp} OUTPUT={outp} VALIDATION_STRINGENCY=SILENT """
        """CREATE_INDEX=true TMP_DIR=/tmp""".format(picard = PICARD_JAR,
        inp = in_file_path, outp = out_file_path))
    run_cmd(command_str, out_log_file_path, err_log_file_path)


def samtools_mapq20(input_file, output_file, file_name):
    """

    :param input_file:
    :param output_file:
    :param file_name:
    :return:
    """
    out_file_path = MAPQ20_OUT_DIR + "12429.samtools.MAPQ20" + ".bam"
    err_log_file_path = LOG_DIR + "12429.samtools.MAPQ20" + ".err.log"
    in_file_path = CLEANSAM_OUT_DIR + "12429.picard.CleanSam.bam"

    command_str = ("""samtools view -bq 20 {inp}""".format(
        inp = in_file_path))
    run_cmd(command_str, out_file_path, err_log_file_path)


def picard_fixmate(input_file, output_file, log_file):
    """

    :param input_file:
    :param output_file:
    :param log_file:
    :return:
    """
    pass


def picard_mark_duplicates(input_file, output_file, log_file):
    """

    :param input_file:
    :param output_file:
    :param log_file:
    :return:
    """
    pass

def picard_collect_duplicates(input_file,
                              output_file, chart_output_file,
                              reference_seq,
                              log_file):
    """

    :param input_file:
    :param output_file:
    :param log_file:
    :return:
    """
    pass


def samtools_mpileup(input_files,
                              output_file,
                              reference_seq,
                              log_file):
    """

    :param input_file:
    :param output_file:
    :param log_file:
    :return:
    """
    pass


def popoolation2_mpileup_to_sync(input_file, output_file, log_file):
    """
    :param input_file:
    :param output_file:
    :param log_file:
    :return:
    """
    mpileup_path = POPOOLATION2_DIR + "mpileup2sync.jar"
    in_file_path = ("/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.mpileup")
    log_file_path = "/share/volatile_scratch/kzukowski/pgi/cattle/test/log/12766.mpileup2sync.log"
    out_file_path = "/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.sync"

    command_str = ("""java -Xmx8g -jar {mpileup} --input {inp}
--output {outp} --fastq-type sanger --min-qual 20 --threads 2"""
        """CREATE_INDEX=true TMP_DIR=/tmp""".format(mpileup = mpileup_path, inp = in_file_path, outp = out_file_path))

    logs = run_cmd(command_str)
    #write to log file
    print(logs)
    pass

def get_all_init_filepaths(dir_path):
    """This function gets the BAM files from the INIT DIR folder.
    """
    init_bam_files = [(dir_path + f) for f in os.listdir(dir_path)
                        if os.path.splitext(f)[1] == '.bam']
    return init_bam_files



if __name__ == '__main__':

    ensure_path_exists(BASE_DIR)
    ensure_path_exists(CLEANSAM_OUT_DIR)
    ensure_path_exists(LOG_DIR)
    print get_all_init_filepaths(INIT_DIR)

    #picard_cleansam(1,2,3)
    samtools_mapq20(1,2,3)
    sys.exit(0)






