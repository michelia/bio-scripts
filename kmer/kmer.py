#encoding=utf8
from __future__ import division
import sys
sys.path.append('/scgene/tiger/invent/guoshuguang/repo_michelia/')
from michelia import SeqIO, dump, load, path

''' 把所有程序写道一个文件中， 相当与总结吧'''

def seq2kmer(seq, kmerLen=17):
    seqLen = len(seq)
    i = 0
    while True:
        if i+kmerLen > seqLen:
            break
        yield seq[i:i+kmerLen]
        i += 1

def fastq2kmer(fqFile):
    for record in SeqIO.parse(fqFile, 'fastq'):
        seq = str(record.seq)
        if 'N' in seq:
            continue
        for kmer in seq2kmer(seq):
            # print kmer
            yield kmer

def kmer2num(kmer, kmerLen=17):
    '''
    kmer is seq base 
    GAAAAAAAAAAAAAGGT  >> 12884901949
    '''
    if len(kmer) != kmerLen:
        raise Exception('kmer length is not 17(kmerLen)')
    base2num = {'A': 0, 'T': 1, 'C': 2, 'G':3}
    kmerNum = 0
    for i, base in enumerate(kmer[::-1]):  # reversal kmer
        kmerNum += 4**i * base2num[base]
    return kmerNum 

def num2piso(num):
    '''
    每个piso中包含一个有5*10**5个元素的字典
    '''
    pisoID = num // 500000
    pisoFile = '/scgene/tiger/invent/guoshuguang/kmerDate/piso5/kmers_%s.piso' % pisoID
    if path(pisoFile).exists():
        kmersDict = load(pisoFile)
        kmersDict[num] = kmersDict.get(num, 0) + 1
        dump(kmersDict, pisoFile)
    else:
        dump({num: 1}, pisoFile)

def sta_piso(pisoDic):
    staDict = {}
    for piso in path(pisoDic).files():
        kmersDict = load(piso)
        for time in kmersDict.iterkeys():
            staDict[time] = staDict.get(time, 0) + 1
    return staDict



# fastqFile = '/scgene/tiger/invent/guoshuguang/kmerDate/test2.fq'
fastqFile = '/scgene/tiger/invent/guoshuguang/kmerDate/rice_genome_1_500m.fq'

for kmer in fastq2kmer(fastqFile):
    num = kmer2num(kmer)
    # print num
    num2piso(num)
