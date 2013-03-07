#! /usr/bin/python

import sys
import re
import commands
import time


def getFiles(dir, subfix): #get files by using ls, return file_list
    cmd = 'ls %s/*%s' % (dir, subfix)
    (status, output) = commands.getstatusoutput(cmd)
    if status == 0:
        file_list = output.split('\n')
      	return file_list
    else:
        sys.stderr.write('[%s] ls get files error: id %d, info: %s\n' %\
                (time.ctime(), status, output))
      	sys.exit(1)


def splitSam(sam_dir, output_dir, sample_name):
    is_head = 0
    sam_list = getFiles(sam_dir, 'sam')

    (status, division_list) =\
            commands.getstatusoutput('ls %s/*sam' % (output_dir))
    if status == 0:
        print >>sys.stderr, "sam files in division directory!"
        sys.exit(1)

    align_result = {}
    for sam in sam_list:
	print >>sys.stderr, "division", sam
	if is_head == 0:
	    common_header = []
	    for line in open(sam):
	        line = line.rstrip()
		if re.search('^@SQ', line):
                    chrom = line.split()[1].split(':')[1]
                    if chrom in align_result:
                        print >> sys.stderr, 'chromosome id error', line
                        sys.exit(1)
                    align_result[chrom] = [line]
                elif re.search('^@', line) and not re.search('^@SQ', line):
                    common_header.append(line)
                else:
		    is_head = 1 #finished header detect
                    break

	    for header in common_header:
                if re.search('^@HD', header):
                    for chrom in align_result:
                        align_result[chrom].insert(0, header)
                else:
                    for chrom in align_result:
                        align_result[chrom].append(header)

        block = 0
        SAM = open(sam, 'r')
        while True:
            line0 = SAM.readline()
            line1 = ""
            line0 = line0.rstrip()
            if len(line0) == 0:
                break
            if not re.search("^@", line0):
                line1 = SAM.readline()
                line1 = line1.rstrip()

                samflag0 = int(line0.split()[1])
                samflag1 = int(line1.split()[1])
                (pair0, chrom0) = getsamflag(samflag0)
                (pair1, chrom1) = getsamflag(samflag1)

                if chrom0 == '':
                    chrom0 = line0.split()[2]
                if chrom1 == '':
                    chrom1 = line1.split()[2]

                if pair0 != "pair" or pair1 != "pair" or chrom0 != chrom1:
                    if chrom0 != "unmap":
                        chrom0 = "unpair"
                    if chrom1 != "unmap":
                        chrom1 = "unpair"

                if chrom0 in align_result:
                    align_result[chrom0].append(line0)
                else:
                    align_result[chrom0] = [line0]

                if chrom1 in align_result:
                    align_result[chrom1].append(line1)
                else:
                    align_result[chrom1] = [line1]

                block += 1

            if block == 1000000:
                output(align_result, output_dir, sample_name)
                align_result = {}
                block = 0;
        SAM.close()

    output(align_result, output_dir, sample_name) #end result



def getsamflag(samflag):
    samflag = int(samflag)

    flags = (1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1)
    pair = 'unpair'
    chrom = ''
    for i in range(0, 11): #unmap
        if flags[i] <= samflag:
            samflag -= flags[i]
            if flags[i] == 8 or flags[i] == 4:
                chrom = 'unmap'
            if flags[i] == 2:
                pair = 'pair'

    if pair == 'unpair' and chrom != 'unmap':
        chrom = 'unpair'

    return (pair, chrom)


def output(align_result, output_dir, sample_name):
    for chrom in align_result:
        align_output = open\
                (output_dir + '/' + sample_name + '_' + chrom + '.sam', 'a')
      	align_output.write("%s\n" % ("\n".join(align_result[chrom])))
        align_output.close()


if len(sys.argv) < 2:
    print "python splitSam.py <sam_dir> <output_dir> <sample_name>"
    sys.exit(0)
else:
    splitSam(sys.argv[1], sys.argv[2], sys.argv[3])
