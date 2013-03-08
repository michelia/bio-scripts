#encoding=utf8
from __future__ import division
import sys
sys.path.append('/scgene/tiger/invent/guoshuguang/repo_michelia/')
from michelia import SeqIO, dump, load, path
from itertools import groupby

''' 把所有程序写道一个文件中， 相当与总结吧'''

splitStep = 500000   # 这一个 每一份所包含的最大的kmer个数。

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

def nums2piso(pisoID, nums):
    '''
    每个piso中包含一个有5*10**5个元素的字典
    '''
    pisoFile = '/scgene/tiger/invent/guoshuguang/kmerDate/piso5/kmers_%s.piso' % pisoID
    # pisoFile = '/scgene/tiger/invent/guoshuguang/kmerDate/testpiso/kmers_%s.piso' % pisoID
    if path(pisoFile).exists():
        kmersDict = load(pisoFile)
    else:
        kmersDict = {}
    for num in nums:
        # print num
        kmersDict[num] = kmersDict.get(num, 0) + 1
    dump(kmersDict, pisoFile)

def sta_piso(pisoDic):
    staDict = {}
    for piso in path(pisoDic).files():
        kmersDict = load(piso)
        for time in kmersDict.iterkeys():
            staDict[time] = staDict.get(time, 0) + 1
    return staDict


def numList2piso(numList):
    numList = sorted(numList)
    for pisoID, nums in groupby(numList, key=lambda x: x//splitStep):
        # print pisoID, nums
        nums2piso(pisoID, nums)

# fastqFile = '/scgene/tiger/invent/guoshuguang/kmerDate/test2.fq'
fastqFile = '/scgene/tiger/invent/guoshuguang/kmerDate/rice_genome_1_500m.fq'

# def get_nums(fastqFile)
numList = []
flag = 0
for kmer in fastq2kmer(fastqFile):
    num = kmer2num(kmer)
    numList.append(num)
    flag += 1
    if flag == 500000:
        numList2piso(numList)
        flag = 0
        numList = []
numList2piso(numList)
