import threading,random,time
from multiprocessing import Process


from multiprocessing import Lock, Process, Queue, current_process, pool
import time
import queue # imported for using queue.Empty exception

# a = "teste"
# def do_job(lista,tasks_to_accomplish, tasks_that_are_done):
#     global a
#     while True:
#         try:
#             '''
#                 try to get task from the queue. get_nowait() function will 
#                 raise queue.Empty exception if the queue is empty. 
#                 queue(False) function would do the same task also.
#             '''
#             task = tasks_to_accomplish.get_nowait()
#         except queue.Empty:

#             break
#         else:
#             '''
#                 if no exception has been raised, add the task completion 
#                 message to task_that_are_done queue
#             '''
#             print(task)
#             tasks_that_are_done.put(task + ' is done by ' + current_process().name)
#             a = task
#             time.sleep(.5)
#     return True


# def main():
#     lista = Queue()
#     number_of_task = 10
#     number_of_processes = 4
#     tasks_to_accomplish = Queue()
#     tasks_that_are_done = Queue()
#     processes = []

#     for i in range(number_of_task):
#         tasks_to_accomplish.put("Task no " + str(i))

#     # creating processes
#     for w in range(number_of_processes):
#         p = Process(target=do_job, args=(lista,tasks_to_accomplish, tasks_that_are_done))
#         processes.append(p)
#         p.start()

#     # completing process
#     for p in processes:
#         p.join()

#     # print the output
#     while not tasks_that_are_done.empty():
#         print(tasks_that_are_done.get())
#     print(a)
#     return True


# if __name__ == '__main__':
#     main()

# n = 0

from multiprocessing import Pool,Queue, Manager
# from queue import Queue

import time

work = (["A", 4], ["B", 2], ["C", 1], ["D", 3])



lista = []

def work_log(work_data):
    # teste = []
    print(" Process %s waiting %s seconds" % (work_data[0], work_data[1]))
    time.sleep(int(work_data[1]))
    print(" Process %s Finished." % work_data[0])
    return work_data[1]
def pool_handler():
    p = Pool(4)
    resultados = p.map(work_log, work)
    # teste = Manager().Queue()
    p.close()
    p.join()
    print(resultados)

if __name__ == '__main__':
    pool_handler()
    # while not teste.empty():
    #     lista.append(teste.get())


    print(lista)

# listaCompleta = []

# def func(lista,numI,numF):
#     for j in range(numI,numF):
#         lista.put(j)

# def main():
#     global listaCompleta

#     lista = Queue()
#     lastNumber = 0
#     procs = []

#     for process in range(1,5):
#             proc = Process(target=func,args=(lista,lastNumber,100*process))
#             proc.start()
#             lastNumber = 100*process
#             procs.append(proc)
#     for process in procs:
#         process.join()
#     while not lista.empty():
#         listaCompleta.append(lista.get())

# if __name__ == "__main__":
#     main()
    

#     print("Testando")
#     print(listaCompleta)

    


# def coder(number):
#     global n
#     a = 0
#     while a < random.randint(1,10000):
#         a += 1
#     dicionario[number] = a
#     # print ('Coders:  %s' % number)
#     n += 1
#     return

# dicionario = {}

# # def a1 ():
# #     return 1
# # def a2 ():
# #     return 2
# # def a3 ():
# #     return 3
# # H = [a1,
# #      a2, 
# #      a3
# #     ]

# threads = {}

# for k in range(5000):
#     t = threading.Thread(target=coder,args=(k,))
#     # threads[k] = 
#     t.start()
# # print(threading.active_count())
# # initial_thread = threading.active_count()
# # while threading.active_count() != 1:
# #     print(threading.active_count())
# #     if threading.active_count() != initial_thread:
# #         initial_thread = threading.active_count()
# #     pass
# print(dicionario)
# print(n)
# print("qualquer ocisa aí")