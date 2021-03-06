#encoding=utf8
from __future__ import division
import sys, pdb, argparse, time, re
b = pdb.set_trace
sys.path.append('/scgene/tiger/invent/guoshuguang/repo_michelia/')
from michelia import CsvReader, path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_fun(dataDir, outPic):
    
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    # colors = ['#DC143C','#8B008B','#00008B','#483D8B','#00FFFF','#008080','#9400D3','#FF00FF','#66CDAA','#008000']
    colors = ['#8B0000', '#FF8C00', '#008000', '#ADFF2F', '#0000CD', '#000000', '#DC143C', '#E9967A', '#808080', '#FF00FF',]
    files = sorted(path(dataDir).files('*.histo'))
    labels = ('1', '1/2', '1/4')
    for i, dataFile in enumerate(files):
        filename = path(dataFile).name
        num = re.findall(r'\d+', filename)[0]
        x, y = get_data(dataFile)
        tx, ty = get_max(x, y)
        ax.plot(tx, ty, marker='s', ms=7)
        print tx, ty
        ax.annotate(str(tx), xy=(tx, ty), 
                xytext=(tx+5, ty+1000000),
                arrowprops=dict(fc='black', shrink=0.15, width=0.7, headwidth=0.9))
        ax.plot(x, y, linestyle='-', marker='.', label=labels[i], c=colors[i])
    ax.legend()
    ax.set_xlabel('Depth')
    ax.set_ylabel('Frequency')
    ax.set_xlim(0, 80)
    ax.set_ylim(0, 3*10**7)
    fig.savefig(outPic)

def get_max(x, y):
    ymax = max(y[3:])
    ymaxX = y.index(ymax) + 1
    return ymaxX, ymax


def get_data(dataFile):
    reader = CsvReader(dataFile).reader
    data = []
    for record in reader():
        data.append(map(int, record[0].split()))
    return zip(*data)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dataDir",
                        help="")
    parser.add_argument("outPic",
                        help="")
    args = parser.parse_args()
    plot_fun(args.dataDir, args.outPic)

if __name__ == '__main__':
    start_time = time.clock()
    main()
    print 
    print '    [INFO] Elapsed time: %.4f' % (time.clock() - start_time)
    print