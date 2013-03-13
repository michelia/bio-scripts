#encoding=utf8
from __future__ import division
import sys, pdb, argparse, time
b = pdb.set_trace
sys.path.append('/scgene/tiger/invent/guoshuguang/repo_michelia/')
from michelia import CsvReader
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_fun(dataFile, outPic):
    reader = CsvReader(dataFile).reader
    data = []
    for record in reader():
        data.append(map(int, record[0].split()))
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    x, y = zip(*data)
    ax.plot(x, y, linestyle='-', marker='.', mfc='r')
    ax.set_xlabel('Depth')
    ax.set_ylabel('Frequency')
    ax.set_xlim(0, 80)
    # ax.set_ylim(0, 50)
    # plt.show()
    fig.savefig(outPic)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dataFile",
                        help="")
    parser.add_argument("outPic",
                        help="")
    args = parser.parse_args()
    plot_fun(args.dataFile, args.outPic)

if __name__ == '__main__':
    start_time = time.clock()
    main()
    print 
    print '    [INFO] Elapsed time: %.4f' % (time.clock() - start_time)
    print