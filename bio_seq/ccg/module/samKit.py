#! /usr/bin/python

import sys
import os
import re

class seqObj:
    def __init__(self, seqLen = 0):
        self.length = seqLen
        self.regStart = []
        self.regEnd = []
        self.region = []

    def do_combine(self, index):
        length = len(self.region)
        if length == 1:
            return True
        else:
            (setFrom, setTo) = self.region[index]
            deletMe = None
            for i in xrange(length):
                if i == index:
                    continue
                if setFrom >= self.region[i][0] and setFrom <= self.region[i][1]:
                    self.region[i] = (self.region[i][0], max(self.region[i][1], setTo))
                    deletMe = index; break
                elif setTo >= self.region[i][0] and setTo <= self.region[i][1]:
                    self.region[i] = (min(self.region[i][0], setFrom), self.region[i][1])
                    deletMe = index; break
                elif setFrom < self.region[i][0] and setTo > self.region[i][1]:
                    self.region[i] = (setFrom, setTo)
                    deletMe = index; break
        return deletMe

    def combine(self, reg):
        (setFrom, setTo) = reg
        if setFrom > setTo:
            setFrom, setTo = setTo, setFrom
            if setTo > self.length:
                self.length = setTo
        self.region.append((setFrom, setTo))
        if len(self.region) == 1:
            return True
        else:
            while True:
                delet = None
                length = len(self.region)
                for i in xrange(length-1, -1, -1):
                    delet = self.do_combine(i)
                    if delet != None:
                        del self.region[delet]
                        break
                if delet == None:
                    break

    def match(self, reg):
        (setFrom, setTo) = reg
        if setFrom > setTo:
            setFrom, setTo = setTo, setFrom
        if len(self.regStart) > 0:
            start = None
            end = None
            for i in xrange(len(self.regStart)):
                if setFrom <= self.regStart[i]:
                    start = i
                    if setTo >= self.regStart[i]:
                        start = i
                        end = i

                        if setTo > self.regEnd[-1]:
                            end = len(self.regEnd) + 1
                            break

                        for j in xrange(i, len(self.regEnd), 1):
                            if setTo < self.regEnd[j]:
                                end = j + 1
                                setTo = self.regEnd[j]
                                break

                    if end == None:
                        end = start
                    break

            if start == None:# append at end
                self.regStart.append(setFrom)
                self.regEnd.append(setTo)
            elif end == start:# insert
                self.regStart = self.regStart[0:start] + [setFrom] + self.regStart[end:]
                self.regEnd = self.regEnd[0:start] + [setTo] + self.regEnd[end:]
            else:# coverage
                self.regStart = self.regStart[0:start] + [setFrom] + self.regStart[end:]
                self.regEnd = self.regEnd[0:start] + [setTo] + self.regEnd[end:]
        else:# Initial
            self.regStart.append(setFrom)
            self.regEnd.append(setTo)
    def stat(self):
        #res = []
        #for i in xrange(len(self.regStart)):
        #    res.append((self.regStart[i], self.regEnd[i]))
        self.region = sorted(self.region)
        print self.region, type(self.region)

def sam2reg(sam_file, reg_times, **argv):
    argv_def = {'avg_ins':0, 'ins_rd':20, 'type':'sampe', 'seq_len':0}
    for i in argv_def:
        if i not in argv:
            argv[i] = argv_def[i]
    if argv['type'] != 'sampe':
        sys.stderr.write("[ERROR] Unsurport: %s\n" % argv['type'])
        sys.exit()

    res = []
    reg_index = 0
    index = 0
    UNIT = 2
    info = range(UNIT)
    avg_ins = 0
    seq = seqObj()
    for line in open(sam_file):
        if line[0] == '@' or line[0] == '#' or line == '':
            continue

        index += 1
        attr = line.split()

        info[(index + 1) % UNIT] = {
                'insert_size':int(attr[8]),
                'match_set':int(attr[3]),
                }

        if index % UNIT != 0:
            continue

        if argv['avg_ins'] != 0:
            min_val = argv['avg_ins'] - argv['ins_rd']
            max_val = argv['avg_ins'] + argv['ins_rd']
            insert_size = abs(info[0]['insert_size'])
            if insert_size < min_val or insert_size > max_val:
                continue

        avg_ins += abs(info[0]['insert_size'])
        reg = (info[0]['match_set'], info[1]['match_set'])
        seq.combine(reg)
        res.append(reg)
        reg_index += 1

        if reg_index == reg_times:
            break
    seq.stat()

    if avg_ins != 0:
        sys.stderr.write("[INFO] Average insersize: %f\n" % (avg_ins/float(reg_index)))
    else:
        sys.stderr.write("[ERROR] Not catch\n")

    return res

if __name__ == '__main__':
    pass
