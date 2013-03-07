################################################################
# Program: qualityCheck.py
# Function: check the quality of fastq files
#   include A T G C N GC% average quality, mode quality, error rate, Q20, Q30
# Author: Chenjiehu
# Version: 1.1
# Date: Tue Aug  2 13:35:44 CST 2011
# Reversion 1.2
#   Edit: Add Q20, Q30
#   Editor: Chen Jiehu
#   Date: Wed Aug 10 16:14:45 CST 2011
# Reversion 1.3
#   Fix a bug: list_cycle
#   Editor: Chen Jiehu
#   Date: Tue Nov 29 10:27:11 CST 2011
# REversion 1.4
#   Support gzip format
#   Add total information of fastq
#   Editor: Chen Jiehu
#   Date: Fri Nov  9 14:27:57 CST 2012
################################################################

import os
import sys
import gzip
import re

def quaStandar(min_qua, max_qua):
    qua_standar = ""
    qua_standar_diff = 0
    if min_qua == 33 and max_qua == 73:
        qua_standar = "Sanger_Phred+33"
        qua_standar_diff = 33
    elif min_qua == 59 and max_qua == 104:
        qua_standar = "Solexa_Solexa+64"
        qua_standar_diff = 64
    elif min_qua == 64 and max_qua == 104:
        qua_standar = "Illumina1.3+_Phred+64"
        qua_standar_diff = 64
    elif min_qua == 66 and max_qua == 105:
        qua_standar = "Illumina1.5+_Phred+64"
        qua_standar_diff = 64
    elif min_qua == 33 and max_qua == 74:
        qua_standar = "Illumina1.8+_Phred+33"
        qua_standar_diff == 33
    elif max_qua == 73:
        qua_standar = "Edit_Sanger_Phred+33"
        qua_standar_diff == 33
    elif max_qua == 104:
        qua_standar = "Edit_Illumina1.3+_Phred+64"
        qua_standar_diff = 64
    elif max_qua == 74:
        qua_standar = "Edit_Illumina1.8+_Phred+33"
        qua_standar_diff = 33
    elif min_qua == 33:
        qua_standar = "Edit_Phred+33"
        qua_standar_diff = 33
    elif min_qua >= 33 and min_qua < 66:
        qua_standar = "Unknown_Phred+33"
        qua_standar_diff = 33
    else:
        print >>sys.stderr, "QUALITY ERROR, unknown quality standar, raw max\
quality is %d and min quality is %d" % (max_qua, min_qua)
        sys.exit(1)
    return (qua_standar, qua_standar_diff)


def errorRate(quality, qua_diff):
    e_rate = 10 ** (quality / -10.0)
    error_rate = e_rate / (1 + e_rate) * 100
    if qua_diff == 33:
        error_rate = (10 ** (quality / -10.0)) * 100
    return error_rate


def staFQ(fq, out_dir):
    name = os.path.basename(fq)
    dict_qua = {} #dict_qua[cycle][quality]
    dict_gc = {'A' : {}, 'T' : {}, 'G' : {}, 'C' : {}, 'N' : {}}
    list_cycle = []
    all_reads = 0
    raw_min_qua = 100
    raw_max_qua = 0
    if re.search('gz$', name):
        FQ = gzip.open(fq, 'r')
    else:
        FQ = open(fq, 'r')
    while True:
        id_1 = FQ.readline()
	if len(id_1) == 0:
            break
	seq = FQ.readline()
	seq = seq.rstrip()
	id_2 = FQ.readline()
	qua = FQ.readline()
	qua = qua.rstrip()

        all_reads += 1

	length = len(seq)
	if length > len(list_cycle): #initialize cycles list
            #name a list of cycles to store information #Reversion 1.3 edit
            for x in range(len(list_cycle), length):
                list_cycle.append(0)

	for i in range(0, length):
            list_cycle[i] += 1
            base = seq[i]
            base_qua_chr = qua[i]
            base_qua = ord(base_qua_chr)
            if base_qua > raw_max_qua:
                raw_max_qua = base_qua
            if base_qua < raw_min_qua:
                raw_min_qua = base_qua

	    if dict_qua.has_key(i):
                if dict_qua[i].has_key(base_qua):
                    dict_qua[i][base_qua] += 1
                else:
                    dict_qua[i][base_qua] = 1
            else:
                dict_qua[i] = {}

	    if dict_gc[base].has_key(i):
                dict_gc[base][i] += 1
	    else:
                dict_gc[base][i] = 1
    FQ.close()

    (qua_standar, qua_diff) = quaStandar(raw_min_qua, raw_max_qua)

    q20 = 0
    q30 = 0
    all_bases = 0
    all_quality = 0
    all_gc = 0
    all_base_noN = 0
    all_qua = 0
    all_mod_quas = {}

    outplot = open(out_dir + '/' + name + '.plot.dat', 'w')
    outplot.write('#Cycle\tBase\tA\tT\tG\tC\tN\tA%\tT%\tG%\tC%\tN%\tAvgQ\tModeQ\tErrorRate%\n')
    for m in range(0, len(list_cycle)): # do analysis for each cycle
        cycle = m + 1
	bases_per_cycle = float(list_cycle[m])
	if bases_per_cycle == 0:
            print >>sys.stderr, 'error cycle', m
            sys.exit(1)
	all_bases += bases_per_cycle

	A = 0
        T = 0
        G = 0
        C = 0
        N = 0
	#not every cycle has 'ATGCN'
        if dict_gc['A'].has_key(m):
            A = dict_gc['A'][m]
        if dict_gc['T'].has_key(m):
            T = dict_gc['T'][m]
        if dict_gc['G'].has_key(m):
	    G = dict_gc['G'][m]
	if dict_gc['C'].has_key(m):
            C = dict_gc['C'][m]
	if dict_gc['N'].has_key(m):
            N = dict_gc['N'][m]

	A_rate = A / bases_per_cycle * 100
	T_rate = T / bases_per_cycle * 100
	G_rate = G / bases_per_cycle * 100
	C_rate = C / bases_per_cycle * 100
	N_rate = N / bases_per_cycle * 100
        all_base_noN += (A + T + G + C)
        all_gc += (G + C)

	total_qua = 0
	mode_qua = 0
	mode_qua_bases = 0

	for quality in dict_qua[m].keys():
            cquality = quality - qua_diff
            total_qua += (cquality * dict_qua[m][quality])
            if cquality >= 20:
                q20 += dict_qua[m][quality]
            if cquality >= 30:
                q30 += dict_qua[m][quality]
            if  dict_qua[m][quality] > mode_qua_bases:
                mode_qua = cquality
                mode_qua_bases = dict_qua[m][quality]

            if cquality in all_mod_quas:  #for all mod qualitys statics
                all_mod_quas[cquality] += dict_qua[m][quality]
            else:
                all_mod_quas[cquality] = dict_qua[m][quality]

        all_qua += total_qua

	avg_qua = total_qua / bases_per_cycle

        error_rate = errorRate(avg_qua, qua_diff)

	outplot.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\n"\
                % (cycle, bases_per_cycle, A, T, G, C, N, A_rate, T_rate, G_rate, C_rate, N_rate, avg_qua, mode_qua, error_rate))

    q20_rate = q20 / float(all_bases) * 100
    q30_rate = q30 / float(all_bases) * 100
    all_gc_rate = all_gc / float(all_base_noN) * 100
    all_avg_qua = all_qua / float(all_base_noN)
    all_mod_qua = sorted(all_mod_quas.items(), key = lambda d:d[1])[-1][0]
    avg_mod_qua_diff = abs(all_mod_qua - all_avg_qua)
    avg_error_rate = errorRate(all_avg_qua, qua_diff)
    mod_error_rate = errorRate(all_mod_qua, qua_diff)

    outplot.write("\n#File\tTotal_reads\tRead_len\tTotal_base\tGC%\t\
Quality_standar\tQ20\tQ30\tAvg_quality\tMod_quality\tabs(AvgQ-ModQ)\tAvg_error_rate%\t\
Mod_error_rate%\n")
    outplot.write("#%s\t%d\t%d\t%d\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%d\t%.2f\t%.6f\t%.6f\n"\
            % (name, all_reads, len(list_cycle), all_bases, all_gc_rate,\
                qua_standar, q20_rate, q30_rate, all_avg_qua, all_mod_qua,\
                avg_mod_qua_diff, avg_error_rate, mod_error_rate))

    outplot.close()

    plotQua(name, len(list_cycle), out_dir + '/' + name + '.plot.dat', out_dir)
    plotGC(name, len(list_cycle), out_dir + '/' + name + '.plot.dat', out_dir)
    plotER(name, len(list_cycle), out_dir + '/' + name + '.plot.dat', out_dir)

def plotER(baseName, cycles, staticsInfo, out_dir):
    plotConf = open(out_dir + baseName + '.errorRate.conf', 'w')

    xT = '"1" 0.0'
    for i in range(9, cycles, 10):
        xT += ', ' + '"' + str(i+1) + '"' + ' ' + str(i) + '.0'

    conf = '''\
            set terminal svg enhanced size 800 600 fname 'arial' fsize 10
	    set output '%s/%s.errorRate.svg'
	    set boxwidth 0.9 absolute
	    set style fill solid 1.00 border -1
	    set style histogram clustered gap 5 title offset character 0, 0, 0
	    set datafile missing '#'
	    set style data histograms
	    set xtics border in scale 1,0.5 nomirror rotate by 0 offset character 0, 0, 0
	    set xtics (%s)
	    set xtics nomirror
	    set ytics nomirror
	    set yrange [0:0.5]
	    set xlabel "Cycles"
	    set ylabel "Error Rate(X100)"
	    set title "Sequencing Error Rate"
	    plot '%s' using 15 ti "Error Rate"
	    '''\
	    % (out_dir, baseName, xT, staticsInfo)

    plotConf.write(conf)
    plotConf.close()
    cmd = 'gnuplot ' + out_dir + baseName + '.errorRate.conf'
    os.system(cmd)
    #cmd = 'rm ' + out_dir + baseName + '.errorRate.conf'
    #os.system(cmd)

def plotGC(baseName, cycles, staticsInfo, out_dir):
    plotConf = open(out_dir + baseName + '.GC.conf', 'w')

    conf = '''\
            set terminal svg enhanced size 800 600 fname 'arial' fsize 10
	    set output '%s/%s.GC.svg'
	    set style data linespoints
	    set xtics nomirror
	    set ytics nomirror
	    set yrange [0:60]
	    set xrange [0:%d]
	    set xlabel \"Cycles\"
	    set ylabel \"Base Precentage(X100)\"
	    set title \" GC Precentage of Sequencing\"
	    plot '%s' using 1:8 title \"A\", '%s' using 1:9 title \"T\", '%s' using 1:10 title \"G\", '%s' using 1:11 title \"C\", '%s' using 1:12 title \"N\"
	    '''\
	    % (out_dir, baseName, cycles, staticsInfo, staticsInfo, staticsInfo, staticsInfo, staticsInfo)

    plotConf.write(conf)
    plotConf.close()
    cmd = 'gnuplot ' + out_dir + baseName + '.GC.conf'
    os.system(cmd)
    #cmd = 'rm ' + out_dir + baseName + '.GC.conf'
    #os.system(cmd)


def plotQua(baseName, cycles, staticsInfo, out_dir):
    plotConf = open(out_dir + baseName + '.quality.conf', 'w')

    conf = '''\
	    set terminal svg enhanced size 800 600 fname 'arial' fsize 10
	    set output '%s/%s.quality.svg'
	    set style data linespoints
	    set xtics nomirror
	    set ytics nomirror
	    set yrange [0:45]
	    set xrange [0:%d]
	    set xlabel \"Cycles\"
	    set ylabel \"Quality\"
	    set title \"Sequencing Quality\"
	    plot '%s' using 1:13 title \"Average Quality\", '%s' using 1:14 title \"Mode Quality\"
	    '''\
	    % (out_dir, baseName, cycles, staticsInfo, staticsInfo)

    plotConf.write(conf)
    plotConf.close()
    cmd = 'gnuplot ' + out_dir + baseName + '.quality.conf'
    os.system(cmd)
    #cmd = 'rm ' + out_dir + baseName + '.quality.conf'
    #os.system(cmd)

if len(sys.argv) > 1:
    staFQ(sys.argv[1], sys.argv[2])
else:
    print "python qcheck.py <fq> <out_dir>\n"
