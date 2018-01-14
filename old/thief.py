#!/usr/local/bin/python3


import multiprocessing
from multiprocessing import Process, Queue, Pool
import sys

def work(value):
    calc = int(value)**int(value)**5
    calc = str(calc)[0:10]
    #print ('Power 10 of %s is %s' % (value,calc))
    return value,calc

def feeder(workerQueue,alist):
    for i in alist:
        #print(i)
        workerQueue.put(i)
        
def worker(workerQueue,writerQueue):
    
    #if workerQueue.get(block = False) != None:
    #    print ('Screw me')
    #    #for i in range (len(workerQueue.get())):
    #    #    print ('screw me')
    #    #    #par = workerQueue.get(block = False)
    #    #    #print ('Echo from Worker: \n Dealing with:  %s'% (par))
    #
    #else: pass
    #    
    
    #par = workerQueue.get(block = False)
    #while par != None:
    #    try:
    #        print ('Echo from Worker: \n Dealing with: %s'% (par))
    #    except:
    #        break
    #    
    
    while True:
        try:
            par = workerQueue.get(block = False)
            print ('Echo from Worker-- Dealing with: %s'% (par))
            res = work(par)
            print ('**Result is: %s' % (res))
            writerQueue.put((par,res))
            #print(writerQueue.qsize())
        
        except:
            print ('List exhausted, breaking')
            break
    
    
    #while True:
    #    try:
    #        par = workerQueue.get(block = False)
    #        print ('Echo from Worker %s: \n Dealing with:'% (worker,par))
    ##        res = work(par)
    ##        writerQueue.put((par,res))
    #    except:
    #        break
    #print ('**Worker finished **')
#    
def writer(writerQueue, filename):
    print ('\nEntered the writer territory with a queue of length: %s' % (writerQueue.qsize()))
    fh_out = open(filename, "w")
    while True:
        try:
            par, res = writerQueue.get(block = False)
            print ('\nFor the value: %s the power to ten is %s' % (par,res))
            fh_out.write('%s,%s\n' % (par,res))
    #        #print >>fhandle, par, res
    #        print('Echo from Writer - List exhausted, breaking')
        except:
            break
    fh_out.close()
    
if __name__ == '__main__':
    resultfile = './test_result'
    nproc =128
    #nproc = int(multiprocessing.cpu_count()*0.5)
    alist = [1,2,3,4,5,6,7,8,9]
    
    
    ##Method 1 - Sophisticated
    #workerQueue = Queue()
    #writerQueue = Queue()
    #
    #FeedProc = Process(target=feeder, args=(workerQueue,alist))
    #WorkProc = Process(target=worker, args=(workerQueue, writerQueue))
    #WritProc = Process(target = writer, args = (writerQueue, resultfile))
    ##
    #FeedProc.start()
    #FeedProc.join()
    #print ('FeedProc join ends')
    #
    #WorkProc.start()
    #WorkProc.join()
    ##for p in range (nproc):
    ##    #worker = p
    ##    WorkProc.start()
    ##for p in range (nproc):
    ##    WorkProc.join()
    #WritProc.start()
    #WritProc.join()
    
    ###Method 2 ## Takes a list returns a list of results
    npool = Pool(int(nproc))    
    res = npool.map_async(work, alist)
    print(res.get())
    print('These are results: %s' % (res))
    sys.exit()
    

    
  
    