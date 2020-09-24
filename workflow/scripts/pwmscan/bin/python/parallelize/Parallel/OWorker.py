import WorkerInterface as WI
import sys


class OWorker(WI.WorkerInterface):
    """
    A subclass of multiprocessing.Process which can be used as a worker 
    process returning data to the caller process.
    
    About the target function and the Queue
    The instance runs a given target function and is connected to 
    the calling process (such as the main process) through a  multiprocess.Queue.
    This Queue object allows the instance to share results returned by the the target 
    function to the calling process. 
    The OWorker will run the target function once with the arguments given at 
    construction time.

    About synchronization
    A multiprocessing.Lock object is given to the OWorker, at construction. This Lock is 
    used to synchronize logging output between processes.
    Results are sent through a multiprocess.Queue which is process safe.

    Running the target function
    Any function can be given when creating an OWorker instance.
    The target function will be run using the tuple containing the arguments given at construction 
    time. This tuple should contain all the necessary arguments to run the function.
    The results returned by the target function are then put into the queue such that the caller 
    process can access them.
    """
    def __init__(self, queue_out=None, target=None, args=tuple(), lock=None, signal_end=False, debug=False):
        """
        Constructor.
        :param queue_out:  multiprocess.Queue object linking the instance to the calling 
        process. The target function results will be shared with the calling process through 
        this queue.
        :param target: a target function. 
        :param args: a tuple of argument to run the target function.
        :param signal_end: a bool, whether the OWorker is should put a poison pill value in the queue_out
        once it is done working.
        :param lock: a multiprocess.Lock object to allow synchronization of logging messages
        among processes.
        """
        super(OWorker, self).__init__(target=target, lock=lock, debug=debug)
        # a flag, whether a signal should be put in the queue_out to signal the
        # end of work
        self._signal_end = signal_end
        # a queue for the instance to get data from the outside
        self._queue_out = queue_out
        # arguments for the target function
        self._args = args

    # target related methods
    def _get_target(self):
        """
        Overrides WorkerInterface._get_target().
        Gets the target function.
        :return: the target function.
        """
        return super(OWorker, self)._get_target()

    def _get_args(self):
        """
        Gets the tuple containing the target function arguments.
        :return: the argument tuple.
        """
        return self._args

    # queue in related methods

    def _get_queue_out(self):
        """
        Gets the queue.
        :return: the queue.
        """
        return self._queue_out

    # flag related methods

    def _get_signal_end(self):
        return self._signal_end

    def _put_to_queue_out(self, data, timeout=0):
        """
        Puts data in the queue.
        :param data data to put in the queue.
        :param timeout the maximum time in seconds spent trying to put data 
        in the queue before raising a Full exception.
        :return: nothing.
        """
        queue = self._get_queue_out()
        self._write_log(sys.stdout, "_put_to_queue_out() START")
        queue.put(data, timeout=timeout)
        self._write_log(sys.stdout, "_put_to_queue_out() END")

    # lock related methods

    def _get_lock(self):
        """
        Overrides WorkerInterface._get_lock().
        Gets the lock.
        :return: the lock.
        """
        return super(OWorker, self)._get_lock()

    # other

    def _end(self):
        """
        A routine to call once the OWorker ends its work. If _signal_end is True, 
        it puts a poison_pill in the queue_out.
        :return: nothing.
        """
        if self._get_signal_end():
            self._put_to_queue_out(OWorker.get_poison_pill())
        else:
            pass

    def run(self):
        """
        Overrides WorkerInterface._get_lock() which overrides multiprocessing.Process.run().
        Runs the target function, with the given arguments and puts the data in the 
        queue. If signal_end was set to True, it also insert a signal_pill in the queue_out 
        once done working.
        :return: nothing. 
        """
        self._write_log(sys.stdout, "run() started")
        # self._put_to_queue_out(self._target(*self._get_args(), lock=self._get_lock()))
        self._put_to_queue_out(self._target(*self._get_args()))
        # signal end if required
        self._end()
        self._write_log(sys.stdout, "run() ended")
        return
