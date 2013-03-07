#! /usr/bin/python

import sys
import os
import re

import fastKit

def align_lis(subject_lis, query_lis):
    sub_seqs = []
    sub_names = []
    que_seqs = []
    que_names = []

    for i in subject_lis:
        (names, seqs) = fastKit.division(i)
        sub_seqs.extend(seqs)
        sub_names.extend(names)

    for i in query_lis:
        (names, seqs) = fastKit.division(i)
        que_seqs.extend(seqs)
        que_names.extend(names)

    for i in xrange(len(que_seqs)):
        for j in xrange(len(sub_seqs)):
            align(que_seqs[j], sub_seqs[i])

def align(subject_seq, query_seq):
    subject_seq = subject_seq.upper()
    query_seq = query_seq.upper()
    print subject_seq, query_seq
