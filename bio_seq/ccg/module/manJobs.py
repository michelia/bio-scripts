#! /usr/bin/python

import sys
import subprocess
import time

def getstatusoutput(cmd, stdinstr = ''):
    p=subprocess.Popen(cmd, shell=True, universal_newlines=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            )
    #p.stdin.write(stdinstr)
    stdoutdata, stderrdata = p.communicate(stdinstr)
    #p.stdin.close()

    return p.returncode, stdoutdata

def qstatJobs(jobsIds, sleep_time = 10):
    while (len(jobsIds) > 0):
        time.sleep(sleep_time)
        runJobs = []
        for ID in jobsIds:
            ID = ID.rstrip()
            cmd = "qstat -j %s" % ID
            (qstatSta, qsubOutput) = getstatusoutput(cmd)
            if qstatSta == 0:# Running.
                runJobs.append(ID)
            elif qstatSta == 1:# Finish
                sys.stderr.write("[INFO] %s FINISH: %s\n" % (time.ctime(), ID))
            else:
                sys.stderr.write("[INFO] %s ERROR: %s\n" % (time.ctime(), ID))
                sys.exit(1)
        jobsIds = runJobs
    return True

def qsubJobs(shell_addr, soure_alloc, sleep_time = 1):
    shell_addr = shell_addr.rstrip()
    cmd = "qsub %s %s" % (soure_alloc, shell_addr)

    (qsubSta, qsubInfo) = getstatusoutput(cmd)
    qsubID = qsubInfo.split()[2]
    time.sleep(sleep_time)

    if qsubSta == 0:
        sys.stderr.write("[INFO] %s qsub SUCCESS: %s" % (time.ctime(), qsubInfo))
        return qsubID
    else:
        sys.stderr.write("[INFO] %s qsub ERROR: %s" % (time.ctime(), qsubInfo))
        sys.exit(1)

def manJobs(shell_addrs, soure_alloc, sleep_time = 10):
    jobsIds = []
    if type(soure_alloc) != type([]):
        soure_alloc = [soure_alloc for x in shell_addrs]

    for i in xrange(len(shell_addrs)):
        jobsIds.append(qsubJobs(shell_addrs[i], soure_alloc[i]))

    time.sleep(sleep_time)
    qstatJobs(jobsIds, sleep_time)

    return True
