def seq2kmer(seq, kmerLen=17):
    seqLen = len(seq)
    i = 0
    while True:
        if i+kmerLen > seqLen:
            break
        yield seq[i:i+kmerLen]
        i += 1

