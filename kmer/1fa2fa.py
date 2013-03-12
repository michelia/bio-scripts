#encoding=utf8
from __future__ import division
import sys, pdb, gzip, time, argparse
from itertools import imap
b = pdb.set_trace
sys.path.append('/scgene/tiger/invent/guoshuguang/repo_michelia/')
from michelia import SeqIO

def main():
    parser = argparse.ArgumentParser(description='把gz压缩的双末端测序fastq文件，转变成一个fasta文件')
    parser.add_argument("fq1",
                        help="gz 压缩的fastq文件")
    parser.add_argument('fa',
                        help='合并后的fasta文件输出位置')
    args = parser.parse_args()
    paire_fq2fa(args.fq1, args.fa)

def paire_fq2fa(fq1, fa):
    fq1 = gzip.open(fq1, 'rb')
    count = SeqIO.convert(fq1, 'fastq', fa, 'fasta')
    print "  Converted %i records" % count

if __name__ == '__main__':
    start_time = time.clock()
    main()
    print 
    print '    [INFO] Elapsed time: %.4f' % (time.clock() - start_time)
    print