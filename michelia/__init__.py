import sys
sys.path.append('/scgene/tiger/invent/guoshuguang/program/moduls')
from path import path
from Bio import SeqIO
from dict_key_of_max_value import dict_key_of_max_value
# from file_filter import file_filter
from two_dimensional_array_transformation import two_dimensional_array_transformation
from encapsulate_pickle import dump, load
from matplotlibMethod import hzebra, vzebra, adjust_spines, adjust_ticks_label_size
from frange import frange, fxrange
from collections import namedtuple
import michelia_csv
import csv
def csvreader(fileHandel):
    return csv.reader(fileHandel, 'excel-tab')
def csvwriter(fileHandel):
    return csv.writer(fileHandel, 'excel-tab')

# def parse_blastout(blastFilePath):
#     BlastRecorder = namedtuple('BlastRecorder', 'a b')
#     with open(blastFilePath) as blastFile:
#         readers = csvreader(blastFile)
#         for row in readers:
#             yield 
