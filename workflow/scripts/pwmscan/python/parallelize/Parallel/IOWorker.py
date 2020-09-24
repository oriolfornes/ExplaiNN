import IWorker as Iw
import OWorker as Ow
import sys


class IOWorker(Iw.IWorker, Ow.OWorker):

    """ 
    A subclass of multiprocessing.Process which can be used as a worker 
    process getting from the caller process and sending data to the caller 
    process.    
    About the target function and the Queues
    The instance runs a given target function and is connected to the calling process 
    (such as the main process) through two  multiprocess.Queue : a queue_in and a 
    queue_out through which the instance can fetch and send data from and to the caller 
    process.
    This queue_in object allows the instance to get arguments to run the target function. 
    The IOWorker will run the target function as long as the queue_in is not empty and will 
    terminate as soon as all data have been consumed. The results returned by the target 
    function are then put to the queue_out such that the caller process can access them.
    
    About synchronization
    A multiprocessing.Lock object is given to the IOWorker, at construction. This Lock is 
    used to synchronize logging output between processes.
    Argument and results are fetch and send through multiprocess.Queue which are process safe.  
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

    def __init__(self, queue_in=None, queue_out=None, target=None, lock=None, signal_end=False, debug=False):
        """
        Constructor
        :param queue_in: a multiprocess.Queue object linking the instance to the calling 
        process. The target function will be run with arguments taken from this queue and 
        as long as the queue is not empty (queue.empty()). Data in this queue are expected 
        to be tuples containing each time a complete set of arguments for the target function.
        :param queue_out: multiprocess.Queue object linking the instance to the calling 
        process. The target function results will be shared with the calling process through 
        this queue. 
        :param target: a target function to run. It should possess an argument 'lock'.
        :param lock: a multiprocess.Lock object to allow synchronization of the target function.
        """
        # todo solve the diamond inheritance problem, is there a way to declare abstract class inheritance as in C++?
        # because IWorker.__init__() and OWorker.__init__() both call
        # WorkerInterface.__init__(), _debug, _target and _lock are initialized twice.
        # However, it does not matter as long as the same value is passed to both
        # constructors.
        # OWorker._args is left as None
        Iw.IWorker.__init__(self, queue_in=queue_in, target=target, lock=lock, debug=debug)
        Ow.OWorker.__init__(self, queue_out=queue_out, target=target, signal_end=signal_end, lock=lock, debug=debug)

    def run(self):
        """
        Overrides WorkerInterface.run() which overrides multiprocessing.Process.run().
        Runs the target function as long as it can successfully fetch arguments 
        from the queue. If signal_end was set to True, it also insert a signal_pill in 
        the queue_out once done working.
        :return: nothing.
        """
        self._write_log(sys.stdout, "run() started")

        target = self._get_target()
        args = self._get_from_queue_in(timeout=None)

        while not self._is_poison_pill(args):
            # could read from queue_in because it was empty
            if args is None:
                continue
            else:
                self._write_log(sys.stdout, "run() calls target")
                self._put_to_queue_out(target(*args), timeout=None)
                args = self._get_from_queue_in(timeout=None)
        # signal end if required
        self._end()
        self._write_log(sys.stdout, "run() ended")
        return
