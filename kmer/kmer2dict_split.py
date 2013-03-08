#encoding=utf8
from __future__ import division
import sys
sys.path.append('/scgene/tiger/invent/guoshuguang/repo_michelia/')
from michelia import dump, load, path

# kmers0file = '/scgene/tiger/invent/guoshuguang/kmerDate/kmers_0.piso'
# aDict = {}
# for i in xrange(0, 5000000):
#     aDict[i] = 0
# dump(aDict, kmers0file)

# a = load(kmers0file)
# print len(a)

def num2piso(num):
    '''
    每个piso中包含一个有5000000个元素的字典
    '''
    pisoID = num // 5000000
    pisoFile = '/scgene/tiger/invent/guoshuguang/kmerDate/kmer_piso/kmers_%s.piso' % pisoID
    if path(pisoFile).exists():
        kmersDict = load(pisoFile)
        kmersDict[num] = kmersDict.get(num, 0) + 1
        dump(kmersDict, pisoFile)
    else:
        dump({num: 1}, pisoFile)

# num2piso(22)