from seq2kmer import fastq2kmer

from kmer_num import kmer2num

from kmer2dict_split import num2piso

# fastqFile = '/scgene/tiger/invent/guoshuguang/kmerDate/test2.fq'
fastqFile = '/scgene/tiger/invent/guoshuguang/kmerDate/rice_genome_1_500m.fq'

for kmer in fastq2kmer(fastqFile):
    num = kmer2num(kmer)
    # print num
    num2piso(num)

