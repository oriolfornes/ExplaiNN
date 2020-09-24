import multiprocessing as mp
import sys
import abc


class WorkerInterface(mp.Process):
    """
    An abstract class for Worker (multriprocessing.Process subclasses) to implement.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, target=None, lock=None, debug=False):
        super(WorkerInterface, self).__init__()
        # a target function to execute
        self._target = target
        # a lock to give to the target function
        self._lock = lock
        # a flag for debug printing
        self._debug = debug

    # debugging methods

    def _get_debug(self):
        """
        Gets the debug flag.
        :return: the debug flag.
        """
        return self._debug

    def _write_log(self, io_obj, message):
        """
        Write a log message to a given io object if debug is on.
        :param io_obj: an object to write to (a file, stdout, stderr, etc).
        :param message: a message to write.
        :return: nothing
        """
        if self._get_debug():
            self._get_lock().acquire()
            io_obj.write("[Process %s] %s\n" % (str(self.pid), message))
            io_obj.flush()
            self._get_lock().release()

    # overriden methods
    def start(self):
        super(WorkerInterface, self).start()
        self._write_log(sys.stdout, "start()")

    def join(self, timeout=None):
        super(WorkerInterface, self).join(timeout=timeout)
        self._write_log(sys.stdout, "join()")

    # other

    def _is_poison_pill(self, value):
        """
        Checks whether a given value is a poison pill 
        or not.
        :param value: any value to check.
        :return: True of False.
        """
        if value == WorkerInterface.get_poison_pill():
            self._write_log(sys.stdout, "_is_poison_pill() TRUE")
            return True
        else:
            return False

    # abstract methods

    @abc.abstractmethod
    def _get_target(self):
        """
        Returns the target function.
        :return: the target function.
        """
        return self._target

    @abc.abstractmethod
    def _get_lock(self):
        """
        Gets the lock.
        :return: the lock.
        """
        return self._lock

    @abc.abstractmethod
    def run(self):
        """
        Overrides multiprocessing.Process.run().
        Should run the target function.
        :return: 
        """
        pass

    @staticmethod
    def get_poison_pill():
        """
        Returns the poison pill value.
        :return: the poison pill value.
        """
        # I always thought 33 was evil :-)
        return 33