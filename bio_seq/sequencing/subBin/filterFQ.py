##################filterFQ.py################
#filter low quality reads
#Author: Chen Jiehu
#Version: 1.1
#Date: Tue Aug  2 14:03:35 CST 2011
#Reversion: 1.2
#	Edit: Add trim 5 upstream
#	Editor: Chen Jiehu
#	Date: Thu Aug 11 12:12:07 CST 2011
#Reversion: 1.3
#       Edit: Filter low quality algorithm, more strict
#       Editor: Chen Jiehu
#	Date: Tue Jan 17 18:00:09 CST 2012
#Reversion: 1.4
#       Edit: Set min length
#       Editor: Chen Jiehu
#       Date: Mon Jun 18 17:43:51 CST 2012
#Reversion: 1.5
#       Edit: add max low quality base and max N base
#       Editor: Chen Jiehu
#       Date: Fri Feb 22 09:31:33 CST 2013
#############################################

#! /usr/bin/python
import sys
import os
import gzip
import re


def quaConfigue(fq):
    if re.search("gz$", fq):
        FQ = gzip.open(fq, 'r')
    else:
        FQ = open(fq, 'r')

    max_qua = 0
    min_qua = 100

    for i in range(10000):
        rid = FQ.readline()
        if len(rid) == 0:
            break
        seq = FQ.readline()
        symble = FQ.readline()
        quality = FQ.readline()
        quality = quality.rstrip()

        for j in range(len(quality)):
            bqua = ord(quality[j])
            if bqua > max_qua:
                max_qua = bqua
            if bqua < min_qua:
                min_qua = bqua
    if min_qua == 33:
        return 33
    elif min_qua >= 64:
        return 64
    else:
        print >>sys.stderr, "quality check error, min quality is %d" % (min_qua)
        sys.exit(1)


def filterPE(fq_1, fq_2, minQuality, minLength, nnum, lbase, outDir):

    if re.search("gz$", fq_1):
        FQ_1 = gzip.open(fq_1, 'r')
        FQ_2 = gzip.open(fq_2, 'r')
    else:
        FQ_1 = open(fq_1, 'r')
        FQ_2 = open(fq_2, 'r')

    quality_diff = quaConfigue(fq_1)
    print >>sys.stderr, "minmal quality: %s" % (minQuality)
    print >>sys.stderr, "quality difference: +%d" % (quality_diff)

    minQuality = int(minQuality) + quality_diff
    minQuality_chr = chr(minQuality)
    print >>sys.stderr, "min quality ascii: %d, chr: %s" %\
            (minQuality, minQuality_chr)
    minLength = int(minLength)
    print >>sys.stderr, "min length: %d" % (minLength)
    nnum = int(nnum)
    print >>sys.stderr, "max N number: %d" % (nnum)
    lbase = int(lbase)
    print >>sys.stderr, "max low quality base number: %d" % (lbase)


    name_1 = os.path.basename(fq_1)
    name_2 = os.path.basename(fq_2)

    OutFqRemoved1 = gzip.open(outDir + "/" + name_1 + '.removed.gz', 'w')
    OutFqRemoved2 = gzip.open(outDir + "/" + name_2 + '.removed.gz', 'w')
    if re.search("gz$", fq_1):
        OutFqRemained1 = gzip.open(outDir + "/" + name_1, 'w')
        OutFqRemained2 = gzip.open(outDir + "/" + name_2, 'w')
    else:
        OutFqRemained1 = gzip.open(outDir + "/" + name_1 + '.gz', 'w')
        OutFqRemained2 = gzip.open(outDir + "/" + name_2 + '.gz', 'w')

    print >>sys.stderr, "filting %s %s" % (name_1, name_2)

    total_reads_pair = 0
    removed_reads_pair = 0
    remained_reads_pair = 0
    remained_base = 0

    while True:
        id_1 = FQ_1.readline()
        if len(id_1) == 0: #stop while at last line
            break

        total_reads_pair += 1
        seq_1 = FQ_1.readline()
        seq_1 = seq_1.rstrip()
        seq_len = len(seq_1)

        symble_1 = FQ_1.readline()

        quality_1 = FQ_1.readline()
        quality_1 = quality_1.rstrip()

        id_2 = FQ_2.readline()

        seq_2 = FQ_2.readline()
        seq_2 = seq_2.rstrip()

        symble_2 = FQ_2.readline()

        quality_2 = FQ_2.readline()
        quality_2 = quality_2.rstrip()

        if ((id_1.split('/'))[0] != (id_2.split('/'))[0]):
            print >>sys.stderr, "fastq id error !" + id_1 + " -- " + id_2
            sys.exit(1)

	while len(seq_1): #trim fastq 5 upstream low quality reads
	    if seq_1[0] == 'N' or ord(quality_1[0]) <= minQuality:
                seq_1 = seq_1[1:]
                quality_1 = quality_1[1:]
            else:
                break

	while len(seq_2):
            if seq_2[0] == 'N' or ord(quality_2[0]) <= minQuality:
                seq_2 = seq_2[1:]
                quality_2 = quality_2[1:]
            else:
                break

	while len(seq_1): #trim fastq 3 downstream low quality reads
            if seq_1[-1] == 'N' or ord(quality_1[-1]) <= minQuality:
                quality_1 = quality_1[:-1]
                seq_1 = seq_1[:-1]
            else:
                break

        while len(seq_2):
            if seq_2[-1] == 'N' or ord(quality_2[-1]) <= minQuality:
                quality_2 = quality_2[:-1]
                seq_2 = seq_2[:-1]
            else:
                break

	is_remove = 0
        n_count0 = len(seq_1.split("N")) - 1
        n_count1 = len(seq_2.split("N")) - 1
        if n_count0  > nnum or  n_count1 > nnum:
            is_remove = 1

	low_quality_count_1 = 0
	for i in range(len(quality_1)):
	    if ord(quality_1[i]) <= minQuality and seq_1[i] != "N":
                low_quality_count_1 += 1

	low_quality_count_2 = 0
        for j in range(len(quality_2)):
            if ord(quality_2[j]) <= minQuality and seq_2[j] != "N":
                low_quality_count_2 += 1

        if low_quality_count_1 > lbase or low_quality_count_2 > lbase:
            is_remove = 1

	if len(seq_1) < minLength or len(seq_2) < minLength:
            is_remove = 1

	if is_remove == 1:
            removed_reads_pair += 1
            OutFqRemoved1.write(id_1 + seq_1 + "\n+\n" + quality_1 + "\n")
            OutFqRemoved2.write(id_2 + seq_2 + "\n+\n" + quality_2 + "\n")
        else:
            remained_reads_pair += 1
            remained_base += (len(seq_1) + len(seq_2))
            OutFqRemained1.write(id_1 + seq_1 + "\n+\n" + quality_1 + "\n")
            OutFqRemained2.write(id_2 + seq_2 + "\n+\n" + quality_2 + "\n")

    OutFqRemoved1.close()
    OutFqRemoved2.close()
    OutFqRemained1.close()
    OutFqRemained2.close()
    FQ_2.close()
    FQ_1.close()

    remained_rate = remained_reads_pair / float(total_reads_pair)

    print >>sys.stderr, "total reads pair: %d" % (total_reads_pair)
    print >>sys.stderr, "remained reads pair: %d(%.2f)" %\
            (remained_reads_pair, remained_rate)
    print >>sys.stderr, "remained bases: %d" % (remained_base)
    print >>sys.stderr, "removed reads pair: %d" % (removed_reads_pair)

if len(sys.argv) > 1:
	filterPE(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],\
                sys.argv[5], sys.argv[6], sys.argv[7])
else:
	print "python filterFQ.py <fq_1> <fq_2> <minQuality> <minReadLength>\
 <maxNnumber> <maxLowQualityBase> <outDir>"
