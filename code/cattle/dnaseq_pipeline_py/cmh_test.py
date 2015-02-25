__author__ = 'nehiljain'

'''
perl /home/nehil/popoolation2_1201/cmh-test.pl \
--input /share/volatile_scratch/nehil/pgi_wc/mouse/dna_seq/nehil_filtered_sync/filtered_chr_18.sync \
--output /share/volatile_scratch/nehil/pgi_wc/mouse/dna_seq/nehil_cmh_stats_cmh/chr_18.cmh \
--min-count -1 \
--min-coverage -1 \
--max-coverage 1000000000 \
--remove-temp \
--population 3-1,7-2,6-8,5-4 \
>/share/volatile_scratch/nehil/pgi_wc/logs/cmh/filtered_chr_18_sync.log 2>&1

'''


CMH_SCRIPT_PATH = '/home/nehil/popoolation2_1201/cmh-test.pl'
FILTERED_INPUT_DIR = '/share/volatile_scratch/nehil/pgi_wc/mouse/dna_seq/nehil_filtered_sync/'
FILTERED_FILE_PREFIX = 'filtered_chr_'
LOG_DIR = '/share/volatile_scratch/nehil/pgi_wc/logs/cmh/'
CMH_OUTPUT_DIR = '/share/volatile_scratch/nehil/pgi_wc/mouse/dna_seq/nehil_cmh_stats_cmh/'
CMH_OUTPUT_FILE_PREFIX = 'chr_'


def get_all_init_filepaths(dir_path):
    """This function gets the BAM files from the INIT DIR folder.
    """
    pattern = re.compile(".sync");
    init_files = [ f for f in os.listdir(dir_path)
                        if re.search(pattern, f)]
    print "-"*100
    print("INIT Files")
    print(init_files)
    print "-"*100
    return init_files

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


def cmh_test(filename):
    """
    :param input_file: string '12766.sorted.bam'
    :param output_file: string '12766.CleanSam.bam'
    :return: stderror and stdout
    """
    print("filename:", filename)
    out_log_file_path = LOG_DIR + os.path.splitext(filename)[0] + ".out.log"
    err_log_file_path = LOG_DIR + os.path.splitext(filename)[0] + ".err.log"
    out_file_path = output_file_names
    print(filename, out_log_file_path, out_file_path)


    command_str = ("""perl {cmh_pl}  --input {inp}  --output {outp}  """
"""--min-count -1  --min-coverage -1  --max-coverage 1000000000  --remove-temp  """
"""--population 3-1,7-2,6-8,5-4  """
""">{log}""".format(cmh_pl = CMH_SCRIPT_PATH,
        inp = filename, outp = out_file_path, log=out_log_file_path))
    print(command_str)
    # run_cmd(command_str)


def main():
    '''

    :return:
    '''

    init_files = get_all_init_filepaths(FILTERED_INPUT_DIR)
    print("init_files",init_files)
    for fname in init_files:
        cmh_test(fname)
