'make kmer seq'

import itertools
 
kmerfile = 'kmer_seq.txt'
kmerfile = open(kmerfile, 'w')
for kmer in itertools.product('ATCG', repeat=17):
    kmer = ''.join(kmer)
    kmerfile.write('%s\t0\n' % kmer)
kmerfile.close()
