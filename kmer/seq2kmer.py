#encoding=utf8
from __future__ import division
import sys
sys.path.append('/scgene/tiger/invent/guoshuguang/repo_michelia/')
from michelia import SeqIO
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

# a = fastq2kmer('/scgene/tiger/invent/guoshuguang/kmerDate/test2.fq')
# # for kmer in :
# #     print kmer
# print len(list(a))