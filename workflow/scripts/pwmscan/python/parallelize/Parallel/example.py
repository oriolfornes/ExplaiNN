import multiprocessing as mp
from IWorker import IWorker
from OWorker import OWorker
from IOWorker import IOWorker


def task(arg):
    """
    A simple function printing the received argument and returning it.
    :param arg: any number.
    :return: the given argument.
    """
    msg = "Hello from %s" % str(arg)
    print msg
    return arg


if __name__ == "__main__":

    # this lock will be used by the workers function to ensure a synchronized access to
    # stdout (usefull if you use debug=True in workers) -> only one process can print at the time on stdout
    stdout_lock = mp.Lock()

    # IWorker example
    # a queue with the arguments to run task()
    queue_tasks = mp.Queue()
    queue_tasks = mp.Queue()
    for i in xrange(10):
        queue_tasks.put((i,))
    # fill the queue with 1 toxic pill per IWorker, instruct them to terminate
    for i in xrange(4):
        queue_tasks.put(IWorker.get_poison_pill())
    # create the workers which have to run the function task() and take arguments from queue_tasks
    workers = [IWorker(queue_in=queue_tasks, target=task, lock=stdout_lock, debug=False) for _ in xrange(0,4,1)]
    # start all the workers
    for worker in workers:
        worker.start()
    # join the workers to be sure the program waits until they have terminated
    for worker in workers:
        worker.join()
    print "Is queue empty : %s" % str(queue_tasks.empty())
    print "IWorkers done"
    print ""
    print ""


    # OWorker
    queue_results = mp.Queue()
    # create 4 workers which have to run the function task() and take arguments from queue_tasks
    # once done, each worker will put i) the result in the queue and ii) toxic pill indicating it is terminating
    workers = [OWorker(queue_out=queue_results, target=task, args=(10*i,), lock=stdout_lock, signal_end=True, debug=False)
               for i in xrange(0,4,1)]
    # start all the workers
    for worker in workers:
        worker.start()
    # check what is in the result queue
    n_terminated = 0
    while not queue_results.empty() and n_terminated != 4:
        value = queue_results.get(2)
        if value != OWorker.get_poison_pill():
            print value
        else:
            n_terminated += 1
    # join the workers to be sure the program waits until they have terminated
    for worker in workers:
        worker.join()
    print "Is queue empty : %s" % str(queue_tasks.empty())
    print "OWorkers done"
    print ""
    print ""

    # IOWorker
    # create 4 workers which have to run the function task() and take arguments from queue_tasks
    # once done, each worker will put i) the result in the queue and ii) toxic pill indicating it is terminating
    queue_tasks = mp.Queue()
    queue_results = mp.Queue()
    workers = [IOWorker(queue_in=queue_tasks, queue_out=queue_results, target=task, lock=stdout_lock, signal_end=True,
                        debug=False)
               for i in xrange(0, 4, 1)]
    # fill the task queue
    for i in xrange(10):
        queue_tasks.put((100*i,))
    # fill the task queue with 1 toxic pill per IOWorker, instruct them to terminate
    for i in xrange(4):
        queue_tasks.put(IOWorker.get_poison_pill())
    # start all the workers
    for worker in workers:
        worker.start()
    # check what is in the result queue
    n_terminated = 0
    while not queue_results.empty() and n_terminated != 4:
        value = queue_results.get(2)
        if value != OWorker.get_poison_pill():
            print value
        else:
            n_terminated += 1
    # join the workers to be sure the program waits until they have terminated
    for worker in workers:
        worker.join()
    print "Is queue empty : %s" % str(queue_tasks.empty())
    print "IOWorkers done"