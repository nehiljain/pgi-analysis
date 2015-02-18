#__author__ = 'nehiljain'


"""Functions of all the sheel commands in the pipeline which will be used as tasks in Ruffus.
"""

import sys
import subprocess
import os
import ruffus

from ruffus import *


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
    print(os.listdir(dir_path))
    init_bam_files = [ f for f in os.listdir(dir_path)
                        if os.path.splitext(f)[1] == '.bam']
    print(init_bam_files)
    return init_bam_files

init_files = get_all_init_filepaths("/Users/nehiljain/coderepo/pgi-ngs-analysis/code/cattle/data")
print(init_files)

@transform(init_files, suffix(".bam"),
           ".CleanSam.cam")
def picard_cleansam(init_files, output_file_names):
    """

    """
    print(init_files)
    in_file_path = init_files
    print(["picard input file",  init_files])
    out_log_file_path = os.path.splitext(in_file_path)[0] + ".out.log"
    err_log_file_path = os.path.splitext(in_file_path)[0] + ".err.log"
    out_file_path = output_file_names

    command_str = ("""cat {inp} > {outp}""".format(
        inp = in_file_path, outp = out_file_path))
    run_cmd(command_str, out_log_file_path, err_log_file_path)

    print(output_file_names)


@follows(picard_cleansam)
@transform(picard_cleansam, suffix(".CleanSam.cam"),
           ".Task2.dam")
def task2(input_file, output_file_names):
    """

    """
    print(["task2 input file" , input_file])
    in_file_path = input_file
    out_log_file_path = os.path.splitext(in_file_path)[0] + ".out.log"
    err_log_file_path = os.path.splitext(in_file_path)[0] + ".err.log"
    out_file_path = output_file_names
    command_str = ("""cat {inp} > {outp}""".format(
        inp = in_file_path, outp = out_file_path))
    run_cmd(command_str, out_log_file_path, err_log_file_path)
    print(output_file_names)



pipeline_get_task_names()


os.chdir("/Users/nehiljain/coderepo/pgi-ngs-analysis/code/cattle/data")
# pipeline_printout(output_stream = sys.stdout, forcedtorun_tasks = [picard_cleansam, task2], verbose=10,checksum_level=3,verbose_abbreviated_path=1)
pipeline_run(target_tasks = [picard_cleansam,task2], verbose=1,checksum_level=3,verbose_abbreviated_path=1)






