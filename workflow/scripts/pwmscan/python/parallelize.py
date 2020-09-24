import sys
import optparse
import subprocess
import multiprocessing
from Parallel import IOWorker


def check_options(options_obj):
    """
    Checks the presence of some values in the option fetch from the command 
    line. Exits if an expected value cannot be found.
    :param options_obj: the return value of parser.parse_args()
    :return: nothing
    """

    # no option given (all default) -> print to use -h
    if (options_obj.info_file is None) and \
            (options_obj.n_proc is 1) and \
            (options_obj.script is None):
        print "use -h option for help"
        sys.exit(1)

    # at least one option give,
    options_dict = vars(options_obj)
    if options_dict["info_file"] is None:
        sys.stderr.write("Option error! No information file was given (-i)!\n")
        sys.exit(1)
        # raise RuntimeError("Option error! No information file was given (-i)!")
    # this should never be triggered as there is a default value
    elif options_dict["script"] is None:
        sys.stderr.write("Option error! No script was given (-s)!\n")
        sys.exit(1)
        # raise RuntimeError("Option error! No script was given (-s)!")


def fill_queue(file_info, queue, prog, n_process):
    """
    Creates argument tuples to execute worker_task() using a IOWorker instance and put
    the argument tuples in the given queue such that the IOWorkers can access them.
    Specifically, run_worker() requires two arguments : the address of the script to
    execute and a tuple of argument necessary to execute the script. Thus each tuple 
    in the queue will contain i) the script location <script> and ii) a tuple of
    argument created from one line of <file_info> and necessary to execute <script>. 
    For instance, if one tuple in the queue is ("bar/foo.sh", ("-f opt1", "-g opt2"))
    the system call will be "bar/foo.sh -f opt1 -g opt2".
    Eventually, the queue will contain <n_line> tuples and <n_proc> poison pills where 
    <n_line> is the number of lines in <file_info> and <n_proc> the number of processes
    which will be run in parellel.
    
    CAUTION : the number of poison pill should be at least the number of processes run in 
    parallel! If less than this, some processes will run endlessly.
    
    :param file_info: the path to a tsv file containing, on each line, a complete set of 
    argument to run the script.
    :param queue: the queue to fill which will later be used by IOWorkers.
    :param prog: the path to the program which will be executed.
    :param n_process: the number of processes which will be run in parallel.
    :return: nothing.
    """
    # read argument tuples from the given file and store them in the queue
    with open(file_info, "rt") as f:
        for line in f:
            line = line.rstrip()
            value = (prog, tuple(line.split('\t')))
            queue.put(value)

    # put some poison pills in the queue
    for _ in xrange(0, n_process, 1):
        queue.put(IOWorker.IOWorker.get_poison_pill())


def worker_task(script_sh, opt_tuple):
    """
    Execute a system call to run the given script with all the options 
    given in the <opt_tuple> tuple. For instance if <script> is 'bar/foo.sh' 
    and <opt_tuple> is ("-f opt1", "-g opt2", "opt3"), this function will call
    'bar/foo.sh -f opt1 -g opt2 opt3'.
    :param script_sh: the script to execute.
    :param opt_tuple: the options to execute the script.
    :return: a byte string containing the raw results of the system call.
    """
    command = "%s" % script_sh
    for opt in opt_tuple:
        command += " %s" % opt
    try:
        output = subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError as e:
        raise e
    return output


def worker_simulation(queue):

    while not queue.empty():
        arg_tuple = queue.get()
        if arg_tuple == IOWorker.IOWorker.get_poison_pill():
            continue
        else:
            command = arg_tuple[0]
            opt = " ".join(arg_tuple[1])
            print "%s %s" % (command, opt)


def main(info_file, script, n_proc, simulate):
    """
    The main part.
    :param info_file: a string, the path to the file containing 
    the list of arguments to run <script>.
    :param script: a string, the path to the script to execute.
    :param n_proc: the number of parallel processes on which the 
    work should be distributed.
    :param simulate: a boolean, whether or not the work should be 
    simulated, in which case the system calls will be printed 
    on stdout instead of being called.
    :return: 0.
    """
    # set up things for parallel workers
    lock_print = multiprocessing.Lock()
    queue_tasks = multiprocessing.Queue()
    queue_results = multiprocessing.Queue()
    list_processes = [IOWorker.IOWorker(queue_in=queue_tasks, queue_out=queue_results, target=worker_task,
                                        lock=lock_print, signal_end=True, debug=False) for _ in xrange(0, n_proc, 1)]

    # fill the queue with argument tuples and with poison pills for the workers
    fill_queue(info_file, queue_tasks, script, n_proc)

    # only simulate the execution
    if simulate:
        worker_simulation(queue_tasks)
    # run the workers
    else:
        # start workers
        for process in list_processes:
            process.start()

        # get the results from the results queue
        n_process_terminated = 0
        while n_process_terminated < n_proc:
            results = queue_results.get(2)
            if results == IOWorker.IOWorker.get_poison_pill():
                n_process_terminated += 1
                continue
            sys.stdout.write(results)
        # join workers
        for process in list_processes:
            process.join()

    return 0


if __name__ == "__main__":
    # parses options
    usage = "python parallelize.py -i <foo>"
    parser = optparse.OptionParser(version="v1.0",
                                   usage=usage,
                                   description="This program allows to distribute a repetive work, "
                                               "executed by a given script (or program), called  each time "
                                               "with a different set of arguments, on several CPU cores "
                                               "(using separated processes). "
                                               "The script/program arguments should be executable using arguments "
                                               "(e.g. myscript arg1 arg2 or myscript -i arg1 -f arg2). "
                                               "The list of arguments should written in a tab separated file. The "
                                               "file format should be the following : each line is expected "
                                               "to contain a complet list of arguments for a script/program call. "
                                               "The script will be executed as many times as there are lines in "
                                               "the file, and every time with a set of argument corresponding to"
                                               "one line in the file. For instance, if the file contains a line "
                                               "'0.001\\tgreen\\tmyfile.dat' then the following command will be "
                                               "executed 'myscript 0.001 green myfile.dat'. If the file contains "
                                               "'-i 0.001\\t-f green\\tmyfile.dat' then the command will be "
                                               "'myscript -i 0.001 -f green myfile.dat'. "
                                               "This program should NEVER be used to execute program "
                                               "coming from an untrusted source and care should be taken to "
                                               "avoid that two simultaneous/subsequent script/program execution "
                                               "overwrite a same file."
                                               "If the script/program execution returns something through stdout, this "
                                               "will be captured and printed by the main process. The only guarantees "
                                               "regarding the results printing are that the results from each "
                                               "script/program call will be printed at once and that this is done in a "
                                               "synchronized manner (process safe way). "
                                               "If the number of task to execute  is bigger than the number "
                                               "of dedicated CPUs, the tasks are enqueued such that there will "
                                               "never be more CPU working than requested.",
                                   epilog="Written by Romain Groux, May 2017")
    parser.add_option("-i", "--infofile",
                      dest="info_file",
                      type="string",
                      default=None,
                      help="A tsv file containing a list of valid arguments on each row to run the script/program.",
                      metavar="file")
    parser.add_option("-p", "--processes",
                      dest="n_proc",
                      type="int",
                      default=1,
                      help="the number of processes to run in parallel, by default 1.",
                      metavar="n")
    parser.add_option("-s", "--script",
                      dest="script",
                      type="string",
                      default=None,
                      help="the script containing all the processing steps to perform on each individual file.",
                      metavar="script")
    parser.add_option("--simulate",
                      dest="simulate",
                      action="store_true",
                      default=False,
                      help="Simulates the run and only prints the commands which would have been otherwise executed.")

    (options, args) = parser.parse_args()

    # check if the options are good, exit otherwise
    check_options(options)

    # runs the jobs
    code = main(options.info_file, options.script, options.n_proc, options.simulate)
    exit(code)
