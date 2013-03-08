#encoding=utf8
import itertools
import decimal


def frange(start, end, step=1.0):
    step = decimal.Decimal(str(step))
    alist = []
    for i in itertools.count():
        next = start + i * step
        if (step > 0.0 and next >= end) or (step < 0.0 and next <= end):
            break
        alist.append(float(next))
    return alist

def fxrange(start, end, step=1.0):
    step = decimal.Decimal(str(step))
    for i in itertools.count():
        next = start + i * step
        if (step > 0.0 and next >= end) or (step < 0.0 and next <= end):
            break
        yield float(next)

