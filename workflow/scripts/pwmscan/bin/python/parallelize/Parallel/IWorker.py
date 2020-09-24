import WorkerInterface as Wi
import Queue
import sys


class IWorker(Wi.WorkerInterface):
    """
    A subclass of multiprocessing.Process which can be used as a worker 
    process getting data from the caller process.
    
    About the target function and the Queue
    The instance runs a given target function and is connected to 
    the calling process (such as the main process) through a  multiprocess.Queue.
    This Queue object allows the instance to get arguments to run the target function. 
    The IWorker will run the target function as long as the Queue is not empty and will 
    terminate as soon as all data from the Queue have been consumed.
    
    About synchronization
    A multiprocessing.Lock object is given to the IWorker, at construction. This Lock is 
    used to synchronize logging output between processes.
    Argument are fetch through a multiprocess.Queue which is process safe.
    
    Running the target function
    Any function can be given when creating an IWorker instance.
    For each invocation of the target function, the IWorker tries to get a tuple of arguments from 
    the queue_in. This tuple should contain all the necessary arguments to run the function.
    For instance, with a target function f(a,b), the queue should contain (a,b) tuples.
    Once the IWorker starts, the target function is run as long as arguments fetch from the queue 
    are not poison pills (IWorker.get_poison_pill() gives this value and worker._is_poison_pill(value) 
    checks whether a value is a poison pill). Be aware that a IWorker will keep working (alive but not 
    running the target function) even if the queue_in is empty!
    """

    def __init__(self, queue_in=None, target=None, lock=None, debug=False):
        """
        Constructor.
        :param queue_in: a multiprocess.Queue object linking the instance to the calling 
        process. The target function will be run with arguments taken from this queue and 
        as long as the queue is not empty (queue.empty()). Data in this queue are expected 
        to be tuples containing each time a complete set of arguments for the target function .
        :param target: a target function to run.
        :param lock: a multiprocess.Lock object to allow synchronization of logging messages
        among processes.
        """
        super(IWorker, self).__init__(target=target, lock=lock, debug=debug)
        # a queue for the instance to get data from the outside
        self._queue_in = queue_in

    # target related methods
    def _get_target(self):
        """
        Override WorkerInterface._get_target().
        Gets the target function.
        :return: the target function.
        """
        return super(IWorker, self)._get_target()

    # queue in related methods

    def _get_queue_in(self):
        """
        Gets the queue
        :return: the queue
        """
        return self._queue_in

    def _get_from_queue_in(self, timeout=0):
        """
        Try to get data from the queue. If the queue 
        is empty, None is returned.
        :param timeout the maximum time in seconds spent trying to fetch data 
        from the queue before raising an Empty exception.
        :return: data from the queue.
        """
        queue = self._get_queue_in()
        value = None
        try:
            value = queue.get(timeout=timeout)
            self._write_log(sys.stdout, "_get_from_queue_in() SUCCESS")
        # the queue is empty or timeout is reached
        except Queue.Empty:
            self._write_log(sys.stdout, "_get_from_queue_in() FAILURE")
            value = None
        finally:
            return value

    # lock related methods

    def _get_lock(self):
        """
        Override WorkerInterface._get_lock().
        Gets the lock.
        :return: the lock.
        """
        return super(IWorker, self)._get_lock()

    def run(self):
        """
        Overrides WorkerInterface.run() which overrides multiprocessing.Process.run().
        Runs the target function as long as it can successfully fetch arguments 
        from the queue.
        :return: nothing.
        """
        self._write_log(sys.stdout, "run() started")

        args = self._get_from_queue_in(2)
        # work as long as it does not get a poison pill
        while not self._is_poison_pill(args):
            # could read from queue_in because it was empty
            if args is None:
                continue
            else:
                self._write_log(sys.stdout, "run() calls target")
                self._get_target()(*args)
                args = self._get_from_queue_in(2)
        self._write_log(sys.stdout, "run() ended")
        return
