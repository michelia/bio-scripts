#! /usr/bin/python

import sys
import os
import re

import pipeline

REVERSE_MAP = {
        'a':'t', 'A':'T',
        't':'a', 'T':'A',
        'g':'c', 'G':'C',
        'c':'g', 'C':'G',
        'n':'n', 'N':'N',
        }

def reverse(seq):
    res = []
    for i in seq:
        try:
            res.append(REVERSE_MAP[i])
        except:
            res.append(i)
            sys.stderr.write('Unknow base: %s\n' % i)

    #print seq, ''.join(res[::-1]); sys.exit()
    return ''.join(res[::-1])

def openFasta_list(file_lis):
    seqs, names = [], []
    seq, name = [], ''
    for file in file_lis:
        for line in open(file):
            line = line.rstrip()
            if line == '':
                continue
            if line[0] == '>':
                if name != '':
                    seqs.append(''.join(seq))
                    names.append(name)
                name = line
                seq = []
            else:
                seq.append(line)
        if name != '':
            seqs.append(''.join(seq))
            names.append(name)

    return names, seqs

def division(fasta_str):
    lines = fasta_str.split("\n")
    lines = [x.rstrip() for x in lines]
    seqs = []
    names = []
    seq_name = ''
    seq = []
    for line in lines:
        if line == '':
            continue
        if line[0] == '>':
            if seq_name != '':
                names.append(seq_name)
                seqs.append(''.join(seq))
            seq = []
            seq_name = line
        else:
            seq.append(line)
    if seq_name != '':
        names.append(seq_name)
        seqs.append(''.join(seq))

    return names, seqs

def formtFa(fasta_str, **argv):
    def formt(seq):
        seqLen = len(seq)
        seqFmt = []
        for i in xrange(0, seqLen, argv['LEN']):
            if i + argv['LEN'] < seqLen:
                seqFmt.append(seq[i:i+argv['LEN']])
            else:
                seqFmt.append(seq[i:])
        return '\n'.join(seqFmt)

    if not argv.has_key('LEN'):
        argv['LEN'] = 100# Default length
    lines = fasta_str.strip('\n').split('\n')
    res = []
    seq = ''
    seqName = ''

    for line in lines:
        if line[0] == '>':
            if seq != '':
                res.append('%s%s\n' % (seqName, formt(seq)))
            seqName = line.rstrip('\n') + '\n'
            seq = ''
        else:
            seq += line
    if seq != '':
        res.append('%s%s\n' % (seqName, formt(seq)))

    return ''.join(res)

def fqSplit(org_name, org_seq, kmer = 31):
    res = []
    for i in xrange(len(org_seq)):
        seq_name = ">%s#%d" % (org_name, i+1)
        seq = org_seq[i:i+kmer]
        if len(seq) == kmer:
            res.append("%s\n%s\n" % (seq_name, seq))

    return ''.join(res)

class openFq:
    UNITS_LEN = 4
    def __init__(self, fqFile, method = 'r'):
        if method not in ['r', 'w', 'a']:
            sys.stderr.write("[ERROR] Unsupport method: %s\n" % method)
            sys.exit()

        self.method = method
        self.fileName = fqFile
        if re.search(r'\.gz$', fqFile):
            import gzip
            myopen = gzip.open
        else:
            myopen = open
        self.file = myopen(self.fileName, method)
    def __del__(self):
        try:
            self.file.close()
        except:
            pass
    def readline(self):
        res = []
        for i in xrange(self.UNITS_LEN):
            line = self.file.readline().rstrip()
            if line:
                res.append(line)
            else:
                break
        return res
    def next(self):
        res = self.readline()
        if res == []:
            self.close()
            raise StopIteration
        else:
            return res
    def __iter__(self):
        return self
    def close(self):
        self.file.close()

def capFastq(id_file_lis, fastq_file_lis, **argv):
    def readFqId(id_file_lis):
        dic = {}
        for i in id_file_lis:
            for line in open(i):
                try:
                    fqId = argv['REGULER'].match(line).group(1)
                except:
                    continue
                dic[fqId] = 1
        return dic

    if not argv.has_key('REGULER'):# Default reguler for catching the id.
        argv['REGULER'] = re.compile(r'@?([^#\n]+)')
    if not argv.has_key('TYPE'):
        argv['TYPE'] = 'CAP'
    if not argv.has_key('OUTPUT_DIR'):
        argv['OUTPUT_DIR'] = os.getcwd()

    pipeline.makedir(argv['OUTPUT_DIR'])
    FQ_UINT_LIN = 4
    READS = [{'info':range(FQ_UINT_LIN)} for x in range(len(fastq_file_lis))]
    NEWFQ = [{'info':range(FQ_UINT_LIN)} for x in range(len(fastq_file_lis))]
    id_dic = readFqId(id_file_lis)
    id_lis = [x for x in id_dic]
    sys.stderr.write("[INFO] input id num:%d\n" % len(id_lis))

    for i in xrange(len(READS)):# Prepare readin and check output.
        if re.search(r'\.gz$', fastq_file_lis[i]):
            import gzip
            myopen = gzip.open
        else:
            myopen = open
        READS[i]['file'] = myopen(fastq_file_lis[i], 'r')
        NEWFQ[i]['ouput'] = "%s%s%s" % (
                argv['OUTPUT_DIR'].rstrip(os.sep),
                os.sep,
                os.path.basename(fastq_file_lis[i].rstrip('.gz'))
                )
        if os.path.isfile(NEWFQ[i]['ouput']):
            if os.path.samefile(fastq_file_lis[i], NEWFQ[i]['ouput']):
                sys.stderr.write("[ERROR] INPUT can't be OUTPUT: %s\n" %
                        os.path.abspath(NEWFQ[i]['ouput']))
                sys.exit(1)
        NEWFQ[i]['file'] = open(NEWFQ[i]['ouput'], 'w')

    while True:
        for ln_idx in xrange(FQ_UINT_LIN):
            for read in READS:
                read['info'][ln_idx] = read['file'].readline()
        if READS[0]['info'][0] == '':
            break

        fqId = argv['REGULER'].match(READS[0]['info'][0]).group(1)

        if argv['TYPE'] == 'CAP':
            if fqId in id_dic:
                for rd_idx in xrange(len(READS)):
                    NEWFQ[rd_idx]['file'].write(''.join(READS[rd_idx]['info']))
                del id_dic[fqId]

        elif argv['TYPE'] == 'IGNORE':
            if fqId not in id_dic:
                for rd_idx in xrange(len(READS)):
                    NEWFQ[rd_idx]['file'].write(''.join(READS[rd_idx]['info']))

        else:
            sys.stderr.write("[ERROR] not support: %s\n" % argv['TYPE'])
            sys.exit(0)

        if len(id_dic) == 0:
            sys.stdout.write("[INFO] Captrued complite: %d\n" % len(id_lis))
            break

    if len(id_dic) != 0:
        sys.stderr.write("[ERROR] Not find: %d\n" % len(id_dic))

    for i in xrange(len(READS)):# Close file.
        READS[i]['file'].close()
        NEWFQ[i]['file'].close()
        sys.stdout.write("[OUTPUT] %s\n" % os.path.abspath(NEWFQ[i]['ouput']))

def capFastqById(id_lis, fastq_file_lis, ouput_dir = False):
    id_dic = {}
    for i in id_lis:
        id_dic[i] = 1
    P_fqId = re.compile(r'@?([^#\n]+)')

    if len(fastq_file_lis) == 2:
        READ_1 = open(fastq_file_lis[0], 'r')
        READ_2 = open(fastq_file_lis[1], 'r')
        info_1 = range(4)
        info_2 = range(4)
        res = [[], []]
        while True:
            for i in xrange(4):
                info_1[i] = READ_1.readline()
                info_2[i] = READ_2.readline()
            if info_1[0] == '':
                break
            #id_str = info_1[0].rstrip()[1:-2]
            id_str = P_fqId.match(info_1[0]).group(1)
            if id_str in id_dic:
                res[0].append(''.join(info_1))
                res[1].append(''.join(info_2))
                del id_dic[id_str]
                if len(id_dic) == 0:
                    break
        READ_1.close()
        READ_2.close()
        res = [''.join(res[0]), ''.join(res[1])]

        if ouput_dir != False:# Output to file
            ouput_dir = ouput_dir.rstrip('/')
            cap_fq_add = [
                    "%s/%s" % (ouput_dir, fastq_file_lis[0].split('/')[-1]),
                    "%s/%s" % (ouput_dir, fastq_file_lis[1].split('/')[-1]),
                    ]
            OUT_1 = open(cap_fq_add[0], 'w')
            OUT_1.write(res[0])
            OUT_1.close()
            OUT_2 = open(cap_fq_add[1], 'w')
            OUT_2.write(res[0])
            OUT_2.close()

            return cap_fq_add
        else:# Return content, not write file.
            return res
    elif len(fastq_file_lis) == 1:
        pass
    else:
        sys.stderr.write("ERROR: fastq file list:%s", ','.join(fastq_file_lis))

def comSeq(seq1, seq2, seq1_val, seq2_val, rever = False):

    def new_fast(match_index):
        def more_val(val1, val2):
            res = []
            for i in xrange(len(val1)):
                if ord(val1[i]) > ord(val2[i]):
                    res.append(val1[i])
                else:
                    res.append(val2[i])

            return ''.join(res)

        overlap_from = len(seq1) - (match_index + query_len)
        overlap_to = match_index + query_len - 1
        seq = seq1[0:overlap_from] + seq2
        val = seq1_val[0:overlap_from] +\
                more_val(seq1_val[overlap_from:], seq2_val[0:overlap_to + 1]) +\
                seq2_val[overlap_to + 1:]
        return seq, val

    if rever == True:
        seq2 = reverse(seq2)

    query_len = 4
    subject_len = 50
    is_find = False

    query = seq1[-query_len:]
    subject = seq2[0:subject_len]

    match_index = subject.find(query)
    if match_index != -1:
        query2 = seq1[-(query_len + 1):]
        if subject.find(query2) != -1:
            if seq1[len(seq1) - (match_index + query_len):] == \
                    seq2[0:match_index + query_len]:# Must match
                        is_find = True

    if is_find == True:
        overlap = seq2[0:match_index] + query
        (new_seq, new_val) = new_fast(match_index)

        return new_seq, new_val
        #print seq1, seq2, overlap
        #print new_seq, new_val
        #sys.exit()

    return False, False

def seqStat(seqStr, **argv):
    def updataStat(base, baseLen):
        try:
            dic[base + "_MAX_LEN"] = max(dic[base + "_MAX_LEN"], baseLen)
        except:
            dic[base + "_MAX_LEN"] = baseLen

    default = {'NAME':'','FIND':['A','T','G','C','N'], 'STDOUT':True, 'DEBUG':False}
    for i in default:# Initial argv
        if i not in argv:
            argv[i] = default[i]

    seqStr = seqStr.upper()
    length = len(seqStr)
    dic = {}# Inital return dic
    lastBase = ''
    repeatLen = 1

    for i in seqStr:

        if i == lastBase:
            repeatLen += 1
        else:
            updataStat(lastBase, repeatLen)
            lastBase = i
            repeatLen = 1

        try:
            dic[i] += 1
        except:
            dic[i] = 1
    updataStat(lastBase, repeatLen)# Fix last unique base bug

    if argv['DEBUG'] != False:
        for i in dic:
            if i[0] not in argv['FIND']:
                argv['FIND'].append(i[0])

    dic['GC_RATE'] = (dic['G'] + dic['C'])/float(length)
    res = ["%s\tLEN:\t%d\tGC_RATE:\t%f" % (argv['NAME'], length, dic['GC_RATE'])]
    if dic.has_key('N'):
        res[0] += "\tN_TIMES:\t%d" % dic['N']

    for base in argv['FIND']:
        try:
            res.append("%s\tTIMES:\t%d\tMAX_LEN:\t%d" % (
                base,
                dic[base],
                dic[base + "_MAX_LEN"],
                ))
        except:
            pass

    if argv['STDOUT'] != False:
        sys.stdout.write('\n'.join(res) + '\n')

    return dic

def fq2fa(fqInfo, **argv):
    default = {'NAME_SPLIT':'#'}
    for i in default:
        if i not in argv:
            argv[i] = default[i]
    return ">%s%s%s\n%s\n" % (fqInfo[0], argv['NAME_SPLIT'], fqInfo[1], fqInfo[1])

