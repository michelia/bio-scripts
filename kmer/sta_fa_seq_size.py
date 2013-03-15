#encoding=utf8
from __future__ import division
import sys,  pdb, argparse, time
b = pdb.set_trace   #调试
sys.path.append('/scgene/tiger/invent/guoshuguang/repo_michelia/')
from michelia import SeqIO

def main():
    function = ''
    parser = argparse.ArgumentParser(description=function)
    parser.add_argument("fq",
                        help="input fastq")
    args = parser.parse_args()
    
    fq = SeqIO.parse(args.fq, 'fastq')
    readLen = len(fq.next().seq)
    count = 1
    for i in fq:
        count += 1
    print 'fq文件总长度为 %s bp' % (readLen*count)


if __name__ == '__main__':
    start_time = time.clock()
    main()
    print
    print '    [INFO] Elapsed time: %.4f' % (time.clock() - start_time)
    print