#! /usr/bin/python
import sys
import re

def main(up, down, in_fq, out_fq):
	OUT = open(out_fq, 'w')
	start = int(up)
	end = 0 - int(down)
	for line in open(in_fq):
		line = line.rstrip()
		if (not re.search("^@", line)) and len(line) > 1:
			if end == 0:
				line = line[start:]
			else:
				line = line[start : end]
		OUT.write("%s\n" % (line))

	OUT.close()

if len(sys.argv) > 1:
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
else:
	print >>sys.stderr, "python %s <5'length> <3'length> <fq> <filted_fq>" % (sys.argv[0])
	sys.exit(0)
