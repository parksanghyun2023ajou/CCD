# run.py

import os
import signal
import sys
import shlex
import subprocess as sp
import argparse as ap
import multiprocessing as mp
from datetime import datetime
from time import time as timer, sleep

#define paths
dirpos = "../example"
binaryName = "./QCSS_UR"
outpos = "../output"
logpos = "../log"

#
curTime = datetime.now().strftime('%m_%d_%H_%M')
benchlist_dir = []
benchlist     = []
benchDir      = sorted(os.listdir(dirpos))

# result include tuple
benchresult = []

###function define ###
def ExecuteCommand(curBench, result, time_limit):
    inputPath  = "%s/%s" % (dirpos, curBench)
    outputPath = "%s/%s" % (outpos, curBench)
    logPath    = "%s/%s.%s.log" % (logpos, curBench, curTime) 
    exeStr     = "%s %s %s | tee %s" % (binaryName, inputPath, outputPath, logPath)
    print(exeStr)
    run = RunCommand(exeStr, time_limit)
    if(run!=None):
        ourResult = ParseLog(logPath)
        result.append([(curBench[:-5]), ourResult])
        print('Successfully run!')
    else:
        print('Time out!')

def RunCommand( exeStr, time_limit ):
    try:
        proc = sp.Popen(exeStr, shell=True)
        return proc.communicate(timeout = time_limit)
    except sp.TimeoutExpired:
        print('Time out!')
        pid = proc.pid
        proc.terminate()
        proc.kill()
        return None

def ParseLog( logPath ):
    add_2q_num, fidelity, time = 0, 0, 0
    f = open(logPath)
    for line in f:
        if line[0:14] == "# add_2q_num: " :
            add_2q_num = int(line[14:])
        if line[0:12] == "# fidelity: " :
            fidelity = int(line[12:])
        if line[0:11] == "END PROGRAM" :
            time = line.rstrip('\n').split(' ')[-1]
    f.close()
    parse_result = [add_2q_num, fidelity, time]
    return parse_result


def getdirList():
    for curBench in benchDir:
        if curBench[-5:] == ".qasm":
            benchlist_dir.append(curBench)


##############
getdirList()
if len(sys.argv) <= 1 or sys.argv[1] == "help":
    print("usage: python run.py <arg>")
    print("")
    print("list of <arg>")
    print("     help        : print how to use run.py")
    print("     list        : print file number of bench directory")
    print("     all         : execute all bench file in directory")
    print("     #           : execute one bench file(or more) in directory")
    print("                   to look file number do list")
    print("     <benchname> : execute one bench file of benchname")
    sys.exit(1)
elif sys.argv[1] == "list":
    for index, curBench in enumerate(benchlist_dir):
        print( str(index) + " : " + curBench )
    sys.exit(1)
elif sys.argv[1].isdigit() and len(sys.argv) == 2:
    benchNum = int(sys.argv[1])
    try:
        benchlist_dir[benchNum]
    except:
        print("ERROR: NO EXISTING BENCH FILE")
        print("Please do python run.py list")
        sys.exit(1)
    ExecuteCommand( benchlist_dir[benchNum], benchresult, 10000) 
elif sys.argv[1] == "all" or len(sys.argv) >= 3:
        if len(sys.argv) >= 3:
            for argnum in range(1, len(sys.argv)):
                num = int(sys.argv[argnum])
                try:
                    benchlist_dir[num]
                except:
                    print("ERROR: NO EXISTING BENCH FILE")
                    print("Please do python run.py list")
                    sys.exit(1)
                benchlist.append(benchlist_dir[num])
        elif sys.argv[1] == "all":
            benchlist = benchlist_dir

        if __name__ == '__main__':
            procs = []
            manager = mp.Manager()
            results = manager.list()
            for index, curBench in enumerate(benchlist):
                proc = mp.Process(target=ExecuteCommand,
                        args=(curBench, results, 10000)) 
                procs.append(proc)
                proc.start()
            for proc in procs:
                proc.join()
        benchresult = results[:]
else:
    print("ERROR: NO ACCEPTABLE ARGUMENT")
    sys.exit(1)

benchresult = sorted(benchresult, key=lambda row: row[0])

print(" ")
print(" ")
print("----------------------------- Result -------------------------------")
print("                        |                   Our                     |")
print("      bench_name        | add_2q_num  |    fidelity    |   RunTime  |")
print("--------------------------------------------------------------------")
for benchResult in benchresult:
    s1 = '%-23s' % benchResult[0]
    s2 = '%11s' % benchResult[1][0]
    s3 = '%14s' % benchResult[1][1]
    s4 = '%10s' % benchResult[1][2]
    try:
        print( s1, "|", s2, "|", s3, "|", s4, "|") 
    except ZeroDivisionError as e:
        print( s1, e)
print("--------------------------------------------------------------------")