#! /usr/bin/python

import sys
import os
import re

def readsInfo(reads, quality, reads_start, reads_end):
	sub_reads = reads[reads_start - 1 : reads_end]
	
	sub_quality = 0
	for i in range(1, 21):
		base_quality = ord(quality[-i]) - 64
		sub_quality += base_quality

	reads_info = (sub_reads, sub_quality)
	return(reads_info)

def filterDuplication(fq_1, fq_2, output_dir):
	FQ_1 = open(fq_1, 'r')
	FQ_2 = open(fq_2, 'r')

	total_reads = 0
	
	duplication = {}

	while True:
		
		id_1 = FQ_1.readline()
		if len(id_1) == 0: #stop while at last line
			break
		id_2 = FQ_2.readline()
		id_1 = id_1.rstrip()
		id_2 = id_2.rstrip()

		total_reads += 1
	
		seq_1 = FQ_1.readline()
		seq_2 = FQ_2.readline()
		seq_1 = seq_1.rstrip()
		seq_2 = seq_2.rstrip()

		symble_1 = FQ_1.readline()
		symble_2 = FQ_2.readline()

		quality_1 = FQ_1.readline()
		quality_2 = FQ_2.readline()
		quality_1 = quality_1.rstrip()
		quality_2 = quality_2.rstrip()

		(sub_reads_1, sub_qulity_1) = readsInfo(seq_1, quality_1, 30, 70)
		(sub_reads_2, sub_qulity_2) = readsInfo(seq_2, quality_2, 30, 70)

		sub_reads = sub_reads_1 + sub_reads_2
		sub_quality = sub_qulity_1 + sub_qulity_2
		sub_id = id_1.split('/')[0]

		if duplication.has_key(sub_reads):
			duplication[sub_reads].append((sub_id, sub_quality))
		else:
			duplication[sub_reads] = [(sub_id, sub_quality)]
	
	FQ_1.close()
	FQ_2.close()

	unique = {}
	unique_reads = 0

	for key in duplication.keys(): #detected and remove duplication
		dups = len(duplication[key])
		unique_reads += 1
		if dups > 1:

			best_read_id = ''
			max_quality = 0
			for read_id_qua in duplication[key]:
				if max_quality < read_id_qua[1]:
					best_read_id = read_id_qua[0]
					max_quality = read_id_qua[1]

			unique[best_read_id] = 1
			print >> sys.stderr, dups, key, duplication[key]
		else:
			unique[duplication[key][0][0]] = 1
	
	unique_rate = unique_reads / float(total_reads) * 100
	duplication_rate = 100 - unique_rate
	print >> sys.stderr, 'duplication_rate:', duplication_rate
	print >> sys.stderr, 'unique_rate:', unique_rate
	
	FQ_1 = open(fq_1, 'r')
	FQ_2 = open(fq_2, 'r')
	fq_1_name = os.path.basename(fq_1)
	fq_2_name = os.path.basename(fq_2)
	FQ_1_OUTPUT = open(output_dir + '/' + fq_1_name, 'w')
	FQ_2_OUTPUT = open(output_dir + '/' + fq_2_name, 'w')

	while True:
		
		id_1 = FQ_1.readline()
		if len(id_1) == 0:
			break
		id_2 = FQ_2.readline()
		id_1 = id_1.rstrip()
		id_2 = id_2.rstrip()

		total_reads += 1
	
		seq_1 = FQ_1.readline()
		seq_2 = FQ_2.readline()
		seq_1 = seq_1.rstrip()
		seq_2 = seq_2.rstrip()

		symble_1 = FQ_1.readline()
		symble_2 = FQ_2.readline()

		quality_1 = FQ_1.readline()
		quality_2 = FQ_2.readline()
		quality_1 = quality_1.rstrip()
		quality_2 = quality_2.rstrip()
		
		sub_id = id_1.split('/')[0]

		if sub_id in unique:
			FQ_1_OUTPUT.write('%s\n%s\n+\n%s\n' % (id_1, seq_1, quality_1))
			FQ_2_OUTPUT.write('%s\n%s\n+\n%s\n' % (id_2, seq_2, quality_2))

	FQ_1_OUTPUT.close()
	FQ_2_OUTPUT.close()
	FQ_1.close()
	FQ_2.close()
	

if len(sys.argv) < 2:
	print "python duplication.py <fq_1> <fq_2> <output_dir>"
else:	
	filterDuplication(sys.argv[1], sys.argv[2], sys.argv[3])
