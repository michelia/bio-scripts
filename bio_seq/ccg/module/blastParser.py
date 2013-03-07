#! /usr/bin/python

import sys
import re

if __name__ == '__main__':
    sys.stdout.write("[ERROR] This is a parser module for blast\n")
    sys.exit()

P = {}
P['QueryName'] = [re.compile(r'^Query=\s?([^\s]+)'), 1]
P['Length'] = [re.compile(r'^Length=([0-9]+)'), 1]
P['RefName'] = [re.compile(r'^>\s?([^\s]+)'), 1]
P['Query'] = [re.compile(r'^Query\s+([0-9]+)\s+([AaTtGgCcNn-]+)\s+([0-9]+)'), 1, 3, 2]
P['Sbjct'] = [re.compile(r'^Sbjct\s+([0-9]+)\s+([AaTtGgCcNn-]+)\s+([0-9]+)'), 1, 3, 2]
P['Identy'] = [re.compile(r'^\s+Identities = ([0-9]+/[0-9]+) \([0-9]+%\), Gaps = ([0-9]+/[0-9]+) \([0-9]%\)'), 1, 2]

def debug(msg):
    sys.stderr.write(msg + '\n'); sys.exit()

def P_value(re_compile_lis, pro_str):
    if type(re_compile_lis) != type([]):
        re_compile_lis = [re_compile_lis, 1]
    res = []
    for gourp_num in re_compile_lis[1:]:
        try:
            res.append(re_compile_lis[0].search(pro_str).group(gourp_num))
        except:
            pass
    if len(res) == 1:
        res = res[0]
    elif len(res) == 0:
        res = False
    return res

class readBlast:
    def __init__(self, blastFile, method = 'r'):
        if method not in ['r']:
            debug("[ERROR] Unsurported method '%s'" % method)
        self.method = method
        self.fileName = blastFile
        self.file = open(self.fileName, self.method)
        self.nowQueryName = ''
    def close(self):
        self.file.close()
    def __del__(self):
        self.file.close()
    def __iter__(self):
        return self
    def next(self):
        res = self.readinfo()
        if res == False:
            self.close()
            raise StopIteration
        else:
            return res
    def combineBase(self, base_1, base_2):
        dic = {
                'R':['A', 'G'], 'Y':['C', 'T'], 'M':['A', 'C'],
                'K':['G', 'T'], 'S':['G', 'C'], 'W':['A', 'T'],
                }
        if base_1 == 'N' or base_2 == 'N':
            return 'N'
        for i in dic:
            if base_1 in dic[i] and base_2 in dic[i]:
                return i
        debug("[ERROR] combineBase: '%s, %s'" % (base_1, base_2))
    def readinfo(self):
        res = {
                'lines':[],
                'data':{},
                'name':'',
                #'region':[],
                }
        QueryLength = 0
        QueryName = ''
        RefName = ''

        if self.nowQueryName != '':
            QueryName = self.nowQueryName

        # captrue a query record
        for line in self.file:
            if P_value(P['QueryName'], line):
                tmp = P_value(P['QueryName'], line)
                if self.nowQueryName == '':
                    self.nowQueryName = tmp
                    QueryName = tmp
                    continue
                else:
                    self.nowQueryName = tmp# Store the QueryName for next
                    break
            if self.nowQueryName != '':
                res['lines'].append(line)

        if len(res['lines']) == 0:
            return False

        # Splite the query record
        for line in res['lines']:
            if P_value(P['Length'], line) and QueryLength == 0:
                QueryLength = int(P_value(P['Length'], line))
                continue
            if P_value(P['RefName'], line):
                RefName = P_value(P['RefName'], line)
                res['data'][RefName] = {'lines':[]}# Inital
                continue
            if RefName != '':
                res['data'][RefName]['lines'].append(line)# Split
        del res['lines']
        res['name'] = QueryName

        # Parser each record
        for RefName in res['data']:
            res['data'][RefName]['Mismatch'] = {}
            #res['data'][RefName]['Gaps'] = {}
            DO_ANALYSIS = True
            #(Query_from, Query_to, Sbjct_from, Sbjct_to) = 0, 0, 0, 0
            for line in res['data'][RefName]['lines']:
                # Read head information
                if P_value(P['Length'], line):
                    res['data'][RefName]['Length'] = int(P_value(P['Length'], line))
                elif P_value(P['Identy'], line):
                    Identy = P_value(P['Identy'], line)
                    mismatch_rate = Identy[0].split('/')
                    mismatch_rate = float(mismatch_rate[0])/float(mismatch_rate[1])
                    if mismatch_rate < 0.98:###### match well ######
                        DO_ANALYSIS = False
                    else:
                        DO_ANALYSIS = True
                        sys.stderr.write("[MATCH_RATE] %f\n" % mismatch_rate)

                if DO_ANALYSIS == False:
                    continue

                # analysis this record if match well
                elif P_value(P['Query'], line):
                    (self.Query_from, self.Query_to, self.Query_base) = P_value(P['Query'], line)# Query_lis
                    #if Query_from == 0:Query_from = int(self.Query_from)
                elif P_value(P['Sbjct'], line):
                    (self.Sbjct_from, self.Sbjct_to, self.Sbjct_base) = P_value(P['Sbjct'], line)# Sbjct_lis
                    #if Sbjct_from == 0:Sbjct_from = int(self.Sbjct_from)

                    # Analysis mismatch and gaps
                    self.findMismatch(res['data'][RefName]['Mismatch'])
                    #self.findGaps(res['data'][RefName]['Gaps'])

            del res['data'][RefName]['lines']# Free memory
            #Query_to = self.Query_to; Sbjct_to = self.Sbjct_to
            #res['region'].append((
            #    QueryName, "%d_%d" % (Query_from, Query_to),
            #    RefName, "%d_%d" % (Sbjct_from, Sbjct_to)
            #    ))

        return res

    def findMismatch(self, dic):
        self.Query_from = int(self.Query_from)
        self.Query_to = int(self.Query_to)
        self.Sbjct_from = int(self.Sbjct_from)
        self.Sbjct_to = int(self.Sbjct_to)
        if self.Sbjct_from > self.Sbjct_to:
            self.Strand = -1
        else:
            self.Strand = 1
        for i in xrange(len(self.Sbjct_base)):
            if len(self.Sbjct_base) != len(self.Query_base):
                debug("[ERRIOR] %s %s" % (self.Query_base, self.Sbjct_base))
            two_base = [x for x in set([self.Sbjct_base[i].upper(), self.Query_base[i].upper()])]
            if '-' not in two_base and len(two_base) == 2:
                mis = self.Sbjct_from + ((i - self.Sbjct_base[0:i].count('-')) * self.Strand)
                base = self.combineBase(two_base[0], two_base[1])
                try:
                    dic[mis][''] + 1
                except:
                    try:
                        dic[mis][base] = 1
                    except:
                        dic[mis] = {base:1}
                #print dic
                #print mis, self.Sbjct_from, self.Strand, i
                #sys.exit()
        return dic

    def findGaps(self, dic):
        pass
