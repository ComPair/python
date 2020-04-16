#!/usr/bin/env python
"""
------------------------------------------------------------------------
A script that runs the ComPair analysis chain on a multicore system.
Author: Jeremy Perkins (jeremy.s.perkins@nasa.gov)

Modified by Sean Griffin
------------------------------------------------------------------------
"""

import os
def setPath(comPath = ""):
    '''Checks for COMPAIRPATH.  Returns 0 if ok, 1 if bad.'''
    if comPath:
        if '~' in comPath:
            os.environ['COMPAIRPATH'] = os.path.expanduser(comPath)
            return 0
        else:
            os.environ['COMPAIRPATH'] = comPath
            return 0
    elif not ('COMPAIRPATH' in os.environ):
        print('Set or provide COMPAIRPATH')
        return 1
    else:
        return 0

def runCosima(srcFile):

    import subprocess
    import gzip
    import time
    from numpy import random

    lock.acquire()
    
    random_sleep = random.rand() + 1.
    print("Sleeping for {0:f} s. This does nothing if a seed is being used in the simulations.".format(random_sleep))
    time.sleep(random_sleep)
    lock.release()

    print(f"Running cosima on {srcFile}")
    base = os.path.splitext(os.path.basename(srcFile))    
    with open(base[0]+'.cosima.stdout.log', 'wb') as fout:
        with open(base[0]+'.cosima.stderr.log', 'wb') as ferr:
            process = subprocess.call(['cosima', '-s', "120", srcFile],
                                 stdout=fout,
                                 stderr=ferr)

    print(f"Done with {srcFile}!")


def runRevan(simFile, cfgFile, geoFile):

    import subprocess
    import gzip
    print(f"Running revan on {simFile}")

    base = os.path.splitext(os.path.basename(simFile))
    with open(base[0]+'.revan.stdout.log', 'wb') as fout:
        with open(base[0]+'.revan.stderr.log', 'wb') as ferr:
            #p = subprocess.Popen(['revan', '-f', simFile, '-c', cfgFile, '-s', '-a', '-n'],
            p = subprocess.call(['revan', '-a', '-n', '-f', simFile, '-g', geoFile, '-c', cfgFile],
                                stdout=fout,
                                stderr=ferr)

    print(f"Done with {simFile}!")    
    
def runRevan_star(files):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return runRevan(*files)

def getFiles(searchDir = '', extension = 'source'):
    from glob import glob
    if not ('COMPAIRPATH' in os.environ):
        print('Set or provide COMPAIRPATH')
        return ""
    if searchDir:
        return glob(searchDir+'/*.'+extension)
    else:
        return glob(os.environ['COMPAIRPATH']+'/Simulations/PerformancePlotSourceFiles/*.'+extension)

def makeLinks(files, folderName='SimFiles'):
    '''Function to make links in directories.  Probably not useful.'''
    from os import symlink
    from os import chdir
    groups = { 1 : ["Cos0.5","Cos0.7"],
               2 : ["Cos0.6","Cos0.9"],
               3 : ["Cos0.8","Cos1.0"]}
    for filename in files:
        for group,angles in groups.iteritems():
            if any(x in filename for x in angles):
                chdir(folderName + str(group))
                symlink('../'+filename, filename)
                chdir('../')

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
    parser.add_argument("--COMPAIRPATH",help="Path to compair files.  You can set this via an environment variable")
    parser.add_argument("--sourcePath",help="Where the source files live.  If not given, will get from COMPAIRPATH.")
    parser.add_argument("--simPath", help="Where the sim files live (from cosima).")
    parser.add_argument("--revanCfg", help="Revan config file (need full path).")
    parser.add_argument("--geoFile", help="Revan geometry file (need full path).")
    parser.add_argument("--reRun", type=bool, default=False, help="Re run/re write")
    args = parser.parse_args()
    print(args)
    if setPath(args.COMPAIRPATH):
        exit()
    else:
        print(f"COMPAIRPATH set to {os.environ['COMPAIRPATH']}")
    
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
            pool.map(runCosima, srcFiles)

    
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