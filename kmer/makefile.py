# aDict = {}
file1 = open('kmer.data', 'w')
for i in xrange(420000000):
    file1.write('%s\t0\n' % i)

file1.close()
