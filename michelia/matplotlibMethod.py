#encoding=utf8
def vzebra(plt, color='#DCDCDC', alpha=0.5, align='edge', ax=None):
    """
    vertical zebra
    """
    if not ax:
        ax = plt.gca()
    xticklocs = ax.xaxis.get_ticklocs()
    xbar = xticklocs[1::2]
    # ymax = ax.get_ylim()[1] #two methods to get the max 
    ymin, ymax = ax.axis()[2:]
    ax.bar(xbar, [(ymax+abs(ymin))] * len(xbar), width=abs(xticklocs[1] - xticklocs[0]), bottom=ymin,
        color=color, alpha=alpha, align=align, linewidth=0)

def hzebra(plt, color='#DCDCDC', alpha=0.5, align='edge', ax=None):
    """
    horizontal zebra 
    """
    if not ax:
        ax = plt.gca()
    yticklocs = ax.yaxis.get_ticklocs()
    ybar = yticklocs[1::2]
    # xmax = ax.get_xlim()[1] #two methods to get the max 
    xmax = ax.axis()[1]
    ax.barh(ybar, [xmax] * len(ybar), height=abs(yticklocs[1] - yticklocs[0]),
        color=color, alpha=alpha, align=align, linewidth=0)


def adjust_spines(ax, spines, width=None):
    # spines 是保留的spines， width是spines和ticks的大小
    for loc, spine in ax.spines.items():
        if loc in spines:
            # spine.set_position(('outward',10)) # outward by 10 points, 这个就是下图droped的功能， 一般可以注释
            # spine.set_smart_bounds(True)   # 这个是一个附加功能， 一般可以注释掉
            pass
        else:
            spine.set_color('none') # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])
    if width:
        ax.spines['left'].set_linewidth(width)
        ax.spines['right'].set_linewidth(width)
        ax.spines['bottom'].set_linewidth(width)
        ax.spines['top'].set_linewidth(width)
        for line in ax.xaxis.get_ticklines():
            line.set_markersize(width)
            line.set_markeredgewidth(width)
        for line in ax.yaxis.get_ticklines():
            line.set_markersize(width)
            line.set_markeredgewidth(width)

            
def adjust_ticks_label_size(ax, size):
    # 这是调整ticklabels(标签刻度)的大小
    for label in ax.xaxis.get_ticklabels():
        label.set_fontsize(size)
    for label in ax.yaxis.get_ticklabels():
        label.set_fontsize(size)