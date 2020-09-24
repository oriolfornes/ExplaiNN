import os
import sys
import subprocess
import optparse
import tempfile


def check_options(options_obj):
    """
    Checks the arguments passed via the command line.
    Terminates the program if an argument value is not found or is wrong.
    :param options_obj: the return value of parser.parse_args()
    :return: nothing
    """
    # no option given (all default) -> print to use -h
    if (options_obj.file_matrix is None) and \
            (options_obj.seq_files is None) and \
            (options_obj.cutoff is None) and \
            (options_obj.n_proc is 1):
        print "use -h option for help"
        sys.exit(1)

    # at least one option given
    options_dict = vars(options_obj)
    if options_dict["file_matrix"] is None:
        sys.stderr.write("Option error! No matrix file was given (-m)!\n")
        sys.exit(1)
    elif options_dict["seq_files"] is None:
        sys.stderr.write("Option error! No sequence files was given (-s)!\n")
        sys.exit(1)
    elif options_dict["cutoff"] is None:
        sys.stderr.write("Option error! No cutoff score was given (-c)!\n")
        sys.exit(1)
    elif options_dict["n_proc"] < 0 or options_dict["n_proc"] == 0:
        sys.stderr.write("Option error! Incorrect number of processes (-p)!\n")
        sys.exit(1)
    elif (options_dict["windex"] is not None) and (options_dict["windex"] <= 0):
        sys.stderr.write("Options error! Incorrect word index size (-i)!\n")
        sys.exit(1)


if __name__ == "__main__":

    # the dependencies (paths to programs) 
    python_path = "/usr/bin/python" # the python 2.7 interpreter 
    parallelize_path = "/home/local/python/parallelize/parallelize.py" # the parall. dispatcher
    matrix_scan_path = "/home/local/bin/matrix_scan" # matrix_scan
    if not os.path.isfile(python_path):
        sys.stderr.write("Error! Cannot find %s!\n" % python_path)
        sys.exit(1)
    if not os.path.isfile(parallelize_path):
        sys.stderr.write("Error! Cannot find %s!\n" % parallelize_path)
        sys.exit(1)
    if not os.path.isfile(matrix_scan_path):
        sys.stderr.write("Error! Cannot find %s!\n" % matrix_scan_path)
        sys.exit(1)

    # parse options
    parser = optparse.OptionParser(version="v1.0",
                                   description="This program executes matrix_scan in a parallel way. It will "
                                               "spawn one process per sequence file to scan.",
                                   epilog="Written by Romain Groux, October 2017")
    parser.add_option("-m", "--matrix",
                      dest="file_matrix",
                      type="string",
                      help="The file containing the matrix to use.",
                      default=None,
                      metavar="file")
    parser.add_option("-c", "--cutoff",
                      dest="cutoff",
                      type="int",
                      help="The score cutoff.",
                      default=None,
                      metavar="score")
    parser.add_option("-f", "--forward",
                      dest="fwd",
                      action="store_true",
                      help="Scan sequences in forward direction.")
    parser.add_option("-i", "--windex",
                      dest="windex",
                      type="int",
                      help="The word index length, by default 7.",
                      default=None,
                      metavar="widx")
    parser.add_option("-s", "--seqs",
                      dest="seq_files",
                      type="string",
                      default=None,
                      help="A list of sequence files to scan. File should be '\\n' separated (e.g. file1\\nfile2"
                           "\\nfile3). This enables to use direclty the output of a ls command, e.g. -f \"$(ls *fa)\"",
                      metavar="files")
    parser.add_option("-p", "--processes",
                      dest="n_proc",
                      type="int",
                      default=1,
                      help="the number of processes to run in parallel, by default 1.",
                      metavar="n")
    (options, args) = parser.parse_args()
    check_options(options)

    n_proc = options.n_proc
    file_matrix = options.file_matrix
    cutoff = options.cutoff
    windex = options.windex
    fwd = options.fwd
    seq_files = options.seq_files
    seq_files_list = [f for f in seq_files.split('\n')]
    # check the sequence files
    for seq_file in seq_files_list:
        if not os.path.isfile(seq_file):
            sys.stderr.write("Error! Sequence file %s does not exist!\n" % seq_file)
            sys.exit(1)
    # check the matrix file
    if not os.path.isfile(file_matrix):
        sys.stderr.write("Error! Matrix file %s does not exist!\n" % file_matrix)
        sys.exit(1)

    # write a temporary file containing the arguments for each call of the program
    # will be deleted as soon as closed
    file_arg = tempfile.NamedTemporaryFile(mode="w+t")
    windex_str = ""
    if windex is not None:
        windex_str = "-i%s \t" % str(windex)
    fwd_str = ""
    if fwd:
        fwd_str = "-f\t"
    for seq_file in seq_files_list:
        line = "%s%s-m %s\t-c %s\t%s\n" % (fwd_str, windex_str, file_matrix, str(cutoff), seq_file)
        file_arg.write(line)
    file_arg.flush()

    # execute the program by parallel dispatch
    try:
        # for a reason I don't understand, using --simulate option raises "IOError: [Errno 32] Broken pipe"
        command = "%s %s -s %s -i %s -p %s" % \
                  (python_path, parallelize_path, matrix_scan_path, file_arg.name, str(n_proc))
        results = subprocess.check_output(command, shell=True)
        sys.stdout.write(results)
    except subprocess.CalledProcessError as e:
        sys.stderr.write("An error occured during the work dispatch!\n")
        raise e

    sys.exit(0)
