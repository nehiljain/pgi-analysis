#__author__ = 'nehiljain'


"""Functions of all the sheel commands in the pipeline which will be used as tasks in Ruffus.
"""

import sys
import subprocess
import os
import ruffus

PICARD_DIR = "/home/kzukowski/soft/picard-tools-1.119/"
POPOOLATION2_DIR = "/home/kzukowski/soft/popoolation2_1201/"
CATTLE_REF_DIR = "/share/volatile_scratch/kzukowski/pgi/cattle/reference/"
CATTLE_REF_SEQ_FILE = CATTLE_REF_DIR + "Bos_taurus.UMD3.1.dna.toplevel.fa"

def run_cmd(cmd_str):
    """Runs the command as given in command string.

    This function uses subprocess to run shell commands in cmd_str. Throws an exception if run
    command fails.
    returns stdout and stderror in a list"""
    process = subprocess.Popen(cmd_str, stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE, shell = True)
    stdout_str, stder_str = process.communicate()
    if process.returncode != 0:
        raise Exception("Failed to run '%s' \n %s %s Non-zero exit status %s" %
                        (cmd_str, stdout_str, stder_str, process.returncode))
    return [stdout_str, stder_str]


def picard_cleansam(input_file, output_file, log_file):
    """
    :param input_file: string '12766.sorted.bam'
    :param output_file: string '12766.CleanSam.bam'
    :param log_file: string '12766.picard.CleanSam.log'
    :return: stderror and stdout
    """
    cleansam_path = PICARD_DIR + "CleanSam.jar"
    in_file_path = ("""/share/volatile_scratch/nehil/cattle/Bos_taurus_DNAalignment_975/12766/run2332_4/12766.MPS1230"""
        """1871-B04.sorted.bam""")
    log_file_path = "/share/volatile_scratch/kzukowski/pgi/cattle/test/log/12766.picard.CleanSam.log"
    out_file_path = "/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.CleanSam.bam"

    command_str = ("""java -Xmx8g -jar {cleansam} INPUT={inp} OUTPUT={outp} VALIDATION_STRINGENCY=SILENT"""
        """CREATE_INDEX=true TMP_DIR=/tmp""".format(cleansam = cleansam_path, inp = in_file_path, outp = out_file_path))

    logs = run_cmd(command_str)
    #write to log file
    print(logs)



def samtools_mapq20(input_file, output_file, log_file):
    """

    :param input_file:
    :param output_file:
    :param log_file:
    :return:
    """
    pass


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
    in_file_path = ("""//share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.mpileup""")
    log_file_path = "/share/volatile_scratch/kzukowski/pgi/cattle/test/log/12766.mpileup2sync.log"
    out_file_path = "/share/volatile_scratch/kzukowski/pgi/cattle/test/data/12766.sync"

    command_str = ("""java -Xmx8g -jar {mpileup} --input {inp}
--output {outp} --fastq-type sanger --min-qual 20 --threads 2"""
        """CREATE_INDEX=true TMP_DIR=/tmp""".format(mpileup = mpileup_path, inp = in_file_path, outp = out_file_path))

    logs = run_cmd(command_str)
    #write to log file
    print(logs)
    pass




if __name__ == '__main__':
    log = run_cmd("echo 'Test'")
    print(log)
    print(PICARD_DIR, POPOOLATION2_DIR, CATTLE_REF_DIR, CATTLE_REF_SEQ_FILE)
    picard_cleansam(1,2,3)
    sys.exit(0)






