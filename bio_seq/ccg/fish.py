#! /usr/bin/python

import sys

sys.path.append("/scgene/elephant/usr/zhouhaibiao/bin/python/module")
import time
import pipeline
sep = pipeline.os.sep

import sub_functions
F = sub_functions
import sub_hashReads
import sub_contig

USAGE = '''#####################################################
Assemble fish. V1.0.2

USAGE:
    <PROM> <-k kmer> [-i insert_size,sd] [-p 0.9] [-l 100] <-f *.fq> [-r ref.fa] <-o dir>

    options:
    -k    INT       kmer, must be odd.
    -i    INT,INT   library size and insert size standard deviation, two libs input like "200,10 500,15".
    -p    FLOAT     the min probability when a extend can be thinking is right. [0.9]
    -l    INT       sequncing read length.[100]
    -f    FILES     fastq files, order by "-i" parameter input.
    -o    DIR       output direction.

    -r    FILE      seed fasta, program will try to extend every sequnce in this file.

#####################################################'''

def read1fa(faFile):
    seq = []
    for line in open(faFile):
        line = line.rstrip()
        if line[0] == '>':
            name = line
        else:
            seq.append(line)
    return name, ''.join(seq)

def readfas(faFile):
    seq = []
    name = ''
    res = []
    for line in open(faFile):
        line = line.rstrip()
        if line[0] == '>':
            if len(seq):
                res.append((name, ''.join(seq)))
                seq = []
            name = line
        else:
            seq.append(line)
    if len(seq):
        res.append((name, ''.join(seq)))
    return res

def formtFa(seq_str, step = 100):
    res = []
    length = len(seq_str)
    for i in xrange(0, length, step):
        res.append(seq_str[i:i+step])
    return '\n'.join(res)

def main(kmer, insert_size_sd, true_rate, read_len, fqFile, faFile, output_dir):
    argv = {
            'kmer':kmer,
            #'fqFile':[(fqFile_lis, insert_size)],
            #'insert_size':[(insert_size, 10)],
            'fqFile':[],
            'insert_size':[],
            'read_len':read_len,
            'true_rate':true_rate,
            'trans_rate':0.9,
            'librl_rate':0.9,
            'refSeq':[],
            'debug':'debug.txt',
            }

    libs_num = len(insert_size_sd)

    if faFile != False:
        argv['refSeq'] = readfas(faFile)

    if argv['debug']:# debug contigs ouput
        DEBUG = pipeline.open("%s%s%s" % (output_dir, sep, argv['debug']), 'w')

    if libs_num > 0:# parser library detail and push to global var "argv"
        for i in xrange(libs_num):
            (insert_size, sd) = map(int, insert_size_sd[i].split(','))
            argv['insert_size'].append((insert_size, sd))
            argv['fqFile'].append((fqFile[i*2:i*2+2], insert_size))

    OUT_CONTIG = pipeline.open("%s%sk%d.contig" % (output_dir, sep, kmer) , 'w')

    K_graph, K_edges = sub_hashReads.hashReads(argv)# Hashreads to deburjin graph with lib width relationship

    if len(argv['refSeq']) == 0:
        contigs = sub_contig.graph2contig(K_graph, K_edges, argv)# accemble contigs by graph
    else:
        contigs = sub_contig.extendFasta(K_graph, K_edges, argv)# extend contigs by graph

    for name, seq in contigs:# formated result to fasta.
        OUT_CONTIG.write("%s\n%s\n" % (name, formtFa(seq)))

    OUT_CONTIG.close()
    if argv['debug']:
        DEBUG.close()

    return 0

if __name__ == '__main__':
    pars = pipeline.parser(sys.argv[1:], {
        '-k':0,
        '-i':[],
        '-p':0.9,
        '-l':100,
        '-f':[],
        '-r':'NOT_INPUT',
        '-o':'',
        },
        USAGE = USAGE,
        )
    if pars['-r'] == 'NOT_INPUT':
        pars['-r'] = False

    if pars['-k'] > 0:# accept all nums input for testing
        prom_time = time.clock()
        main(pars['-k'], pars['-i'], pars['-p'], pars['-l'], pars['-f'], pars['-r'], pars['-o'])
        prom_time = time.clock() - prom_time
        F.debug("[FINISH] total %.2f time_speed\n" % prom_time)
    else:
        print USAGE
