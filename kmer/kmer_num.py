
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

# print kmer2num('GAAAGAAAAAAAAAGGT')
# num = 12884901949
def num2kmer(num, kmerLen=17):
    '''
    12884901949 >> GAAAAAAAAAAAAAGGT
    '''
    baseNums = (('G', 3), ('C', 2), ('T', 1))
    kmer = ''
    for i in range(kmerLen)[::-1]:
        for base, baseNum in baseNums:
            baseInum =  4**i * baseNum
            if num >= baseInum:
                kmer += base
                num -= baseInum
                break
            elif base == 'T':
                kmer += 'A'
    return kmer
# print num2kmer(num=256)


