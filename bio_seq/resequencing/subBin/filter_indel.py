#! /usr/bin/python
import sys
import re

def main(raw_indel, indel):
    INDEL = open(indel, 'w')
    INDEL.write("CHROM\tPOS\tREF\tALT\tINDEL_QUAL\tReads(most_\
stringent)\tHom/Het\tGQ\tHP\n")
    for line in open(raw_indel):
        line = line.rstrip()
        if re.search("^#", line):
            continue
        (chrom, pos, iid, ref, alt, qua, ffilter, info, fformat, sample) =\
        line.split()
        if ffilter != "PASS":
            continue
        ref_len = len(ref)
        alt_len = len(alt)
        qua = int(qua)
        if ref_len > 10 or alt_len > 10:
            continue
        if qua < 20:
            continue
        (dp, nf, nr, nrs, nfs, hp) = info.split(";")
        dp = int(dp.split("=")[-1])
        nf = int(nf.split("=")[-1])
        nr = int(nr.split("=")[-1])
        nrs = int(nrs.split("=")[-1])
        nfs = int(nfs.split("=")[-1])
        hp = int(hp.split("=")[-1])
        (gt, gq) = sample.split(":")
        gq = int(gq)
        if hp > 10:
            continue
        if nf + nr <= 3:
            continue
        if gq <= 4:
            continue
        hom_het = "HET"
        if gt == "1/1":
            hom_het = "HOM"
        nf_nr = nf + nr
        INDEL.write("%s\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\n" % (chrom, pos, ref,\
                alt, qua, nf_nr, hom_het, gq, hp))
    INDEL.close()


if (len(sys.argv) > 1):
    main(sys.argv[1], sys.argv[2])
else:
    print >>sys.stderr, "python %s <raw_indel> <indel>" % (sys.argv[0])
