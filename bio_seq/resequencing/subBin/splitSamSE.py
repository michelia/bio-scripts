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
                        print >> sys.stderr, 'chromosome id error'
                        sys.exit(1)
                    align_result[chrom] = [line]
		elif re.search('^@', line) and not re.search('^@SQ', line):
                    common_header.append(line)
		else:
                    break

            for header in common_header:
                if re.search('^@HD', header):
                    for chrom in align_result:
                        align_result[chrom].insert(0, header)
                else:
                    for chrom in align_result:
                        align_result[chrom].append(header)
	    is_head = 1 #finished header detect

	block = 0
        for line in open(sam):
            line = line.rstrip()
            block += 1
            if not re.search('^@', line):
                chrom = line.split()[2]
                if chrom == "*":
                    chrom = "unmap"

		if chrom in align_result:
                    align_result[chrom].append(line)
                else:
                    align_result[chrom] = [line]

            if block == 1000000:
                output(align_result, output_dir, sample_name)
                align_result = {}
                block = 0;
	output(align_result, output_dir, sample_name) #end result


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
