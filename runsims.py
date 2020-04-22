#!/usr/bin/env python
"""
------------------------------------------------------------------------
A script that runs the ComPair analysis chain on a multicore system.
Author: Jeremy Perkins (jeremy.s.perkins@nasa.gov)

Modified by Sean Griffin
------------------------------------------------------------------------
"""

import os

def runCosima(srcFile, seed):

    import subprocess
    import gzip
    import time
    from numpy import random

    lock.acquire()
    
    if seed == "":
        random_sleep = random.rand() + 1.
        print("Sleeping for {0:f} s. This is to vary seed times.".format(random_sleep))
        time.sleep(random_sleep)

    lock.release()

    tStart = time.time()
    tshow =time.strftime("%H:%M:%S")

    print(f"Running cosima on {srcFile}; started at {tshow}")
    base = os.path.splitext(os.path.basename(srcFile))    
    with open(base[0]+'.cosima.stdout.log', 'wb') as fout:
        with open(base[0]+'.cosima.stderr.log', 'wb') as ferr:
            if seed != '':
                process = subprocess.call(['cosima', '-s', seed, srcFile],
                                     stdout=fout,
                                     stderr=ferr)
            else: 
                process = subprocess.call(['cosima', srcFile],
                                     stdout=fout,
                                     stderr=ferr)                

    tStop = time.stop()
    
    deltaT = tStop - tStart

    print(f"Done with {srcFile}! Time taken: {deltaT:.1f} seconds. ")


def runRevan(simFile, cfgFile, geoFile):

    import subprocess
    import gzip
    import time 

    
    tStart = time.time()
    tshow =time.strftime("%H:%M:%S")

    print(f"Running revan on {simFile}; started at {tshow}")

    base = os.path.splitext(os.path.basename(simFile))
    with open(base[0]+'.revan.stdout.log', 'wb') as fout:
        with open(base[0]+'.revan.stderr.log', 'wb') as ferr:
            #p = subprocess.Popen(['revan', '-f', simFile, '-c', cfgFile, '-s', '-a', '-n'],
            p = subprocess.call(['revan', '-a', '-n', '-f', simFile, '-g', geoFile, '-c', cfgFile],
                                stdout=fout,
                                stderr=ferr)

    tStop = time.stop()
    deltaT = tStop - tStart
    print(f"Done with {simFile}! Time taken: {deltaT:.1f} seconds.")    
    
def runCosima_star(files):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return runCosima(*files)


def runRevan_star(files):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return runRevan(*files)

def getFiles(searchDir = '', extension = 'source'):
    from glob import glob
    return glob(searchDir+'/*.'+extension)
    
def notDone(sims, tras):
    sims_base = [name[0:-4] for name in sims]
    tras_base = [name[0:-4] for name in tras]
    diff = set(sims_base) - set(tras_base)
    return [name+".sim" for name in diff]

def init_lock(l):
    '''
    Initialize the global multiprocessing lock.
    '''
    global lock
    lock = l

def cli():
    from multiprocessing import Pool, Lock
    from itertools import repeat
    l = Lock()
    helpString = "Submits cosima or revan jobs to multiple cores."
    import argparse
    parser = argparse.ArgumentParser(description=helpString)
    parser.add_argument("jobs", type=int, help="The number of jobs you wish to spawn (usually the number of cores on your machine).")
    parser.add_argument("--runCosima", type=bool, default=False, help="Run cosima (default is false)")
    parser.add_argument("--runRevan", type=bool, default=False, help="Run revan (default is false)")
    parser.add_argument("--sourcePath",help="Where the source files live.  If not given, will get from COMPAIRPATH.")
    parser.add_argument("--simPath", help="Where the sim files live (from cosima).")
    parser.add_argument("--revanCfg", help="Revan config file (need full path).")
    parser.add_argument("--geoFile", help="Revan geometry file (need full path).")
    parser.add_argument("--reRun", type=bool, default=False, help="Re run/re write")
    parser.add_argument("--seed", default='', help="Cosima seed (optional)")    
    args = parser.parse_args()
    print(args)

    if args.runCosima:
        srcFiles = getFiles(args.sourcePath,'source')
        print(srcFiles)
        if not srcFiles:
            print("No source files found")
            exit()
        else:
            print(f"Got this many source files: {str(len(srcFiles))}")
            print("Spawing jobs")
            pool = Pool(processes=args.jobs, initializer=init_lock, initargs=(l,))
            #pool.map(runCosima, srcFiles)
            pool.map(runCosima_star, zip(srcFiles, repeat(args.seed)))

    
    if args.runRevan:
        simFiles = getFiles(args.simPath,'sim')
        for s in simFiles: 
            print(s)
        print("\n\n\n")

        traFiles = getFiles(args.simPath,'tra')
        if not args.reRun:
            simFiles = notDone(simFiles,traFiles)
        if not args.revanCfg:
            print("Need to specify the config file for revan (--revanCfg)")
            exit()
        if not args.geoFile:
            print("Need to specify the geometry file for revan (--geoFile)")
            exit()            
        elif not simFiles:
            print("No sim files found")
            exit()
        else:
            print(f"Got this many sim files: {str(len(simFiles))}")
            print("Spawning jobs")
            pool = Pool(processes=args.jobs)
            pool.map(runRevan_star,
                     zip(simFiles,repeat(args.revanCfg),repeat(args.geoFile)))

if __name__ == '__main__': cli()