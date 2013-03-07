#! /usr/bin/python

import sys
import re
import os

#import other module
sys.path.append("/scgene/elephant/pipeline/python/module")
import manJob
import getFiles
import dirCheck
import makeShell

BIN = sys.argv[0]
PATH = os.path.dirname(BIN)
PATH = PATH.rstrip()
if PATH == "":
    PATH = "."
DINDEL = PATH + "/dindel-1.01-linux-64bit"
MAKEWINDOW = PATH + "/makeWindows.py"
MERGEOPDIPLOID = PATH + "/mergeOutputDiploid.py"
FILTER = PATH + "/filter_indel.py"
SGEHEAD = "#$ -S /bin/sh\n#! /bin/bash"

def stage1(sort_dir, ref_dir, tem_dir):
    stage1_dir = tem_dir + "/stage1"
    dirCheck.mkdir(stage1_dir)
    bam_files = getFiles.get_files(sort_dir, "sort.bam")
    qsub_shells = []
    for bam in bam_files:
        name = os.path.basename(bam)
        base_name = name.split('.')[0]
        chrom = base_name.split('_')[-1]
        shell = "%s\n%s --analysis getCIGARindels --bamFile %s --outputFile %s/%s_dindel_output --ref %s/%s.fa" %\
                (SGEHEAD, DINDEL, bam, stage1_dir, base_name, ref_dir, chrom)
        shell_name = "dindel_%s_stage1.sh" % (base_name)
        qsub_shells.append(makeShell.make_shell(stage1_dir, shell, shell_name))
    manJob.manJobs(qsub_shells, "-cwd -l vf=2G")
    return stage1_dir


def stage2(stage1_dir, tem_dir):
    stage2_dir = tem_dir + "/stage2"
    dirCheck.mkdir(stage2_dir)
    stage1_files = getFiles.get_files(stage1_dir, "variants.txt")
    qsub_shells = []
    for stage1_file in stage1_files:
        name = os.path.basename(stage1_file)
        base_name = "_".join(name.split("_")[:-2])
        shell = "%s\npython %s --inputVarFile %s --windowFilePrefix %s/%s.realign_windows --numWindowsPerFile 1000" %\
                (SGEHEAD, MAKEWINDOW, stage1_file, stage2_dir, base_name)
        shell_name = "dindel_%s_stage2.sh" % (base_name)
        qsub_shells.append(makeShell.make_shell(stage2_dir, shell, shell_name))
    manJob.manJobs(qsub_shells, "-cwd -l vf=2G")
    return stage2_dir


def stage3(sort_dir, ref_dir, stage1_dir, stage2_dir, tem_dir):
    stage3_dir = tem_dir + "/stage3"
    dirCheck.mkdir(stage3_dir)
    bam_files = getFiles.get_files(sort_dir, "sort.bam")
    qsub_shells = []
    for bam in bam_files:
        name = os.path.basename(bam)
        base_name = name.split('.')[0]
        chrom = base_name.split('_')[-1]
        windows = getFiles.get_files(stage2_dir, base_name + ".realign*")
        for window in windows:
            window_name = ".".join(window.split("_")[-1].split(".")[:-1])
            shell = ("%s\n%s --analysis indels --doDiploid --bamFile %s --ref\
 %s/%s.fa --varFile %s --libFile %s/%s_dindel_output.libraries.txt --outputFile\
 %s/%s.%s" % (SGEHEAD, DINDEL, bam, ref_dir, chrom, window,\
                    stage1_dir, base_name, stage3_dir, base_name,\
                    window_name))
            shell_name = "dindel_%s_%s_stage3.sh" % (base_name, window_name)
            qsub_shells.append(makeShell.make_shell(stage3_dir, shell,\
                    shell_name))
    manJob.manJobs(qsub_shells, "-cwd -l vf=1G")
    return stage3_dir


def stage4(ref_dir, stage3_dir, tem_dir):
    stage4_dir = tem_dir + "/stage4"
    dirCheck.mkdir(stage4_dir)
    refs = getFiles.get_files(ref_dir, "chr*.fa")
    qsub_shells = []
    for ref in refs:
        chrom = os.path.basename(ref).split(".")[0]
        windows = getFiles.get_files(stage3_dir, "*" + chrom + ".*.txt")
        sample = "_".join(os.path.basename(windows[0]).split("_")[:-1])
        file_list = stage4_dir + "/" + sample + "_" + chrom + ".list"
        LIST = open(file_list, 'w')
        LIST.write("%s" % "\n".join(windows))
        LIST.close()
        shell = "%s\npython %s --inputFiles %s --outputFile\
 %s/%s_%s.variantCalls.vcf --refFile %s" % (SGEHEAD, MERGEOPDIPLOID,\
                file_list, stage4_dir, sample, chrom, ref)
        shell_name = "dindel_%s_%s_stage4.sh" % (sample, chrom)
        qsub_shells.append(makeShell.make_shell(stage4_dir, shell, shell_name))
    manJob.manJobs(qsub_shells, "-cwd -l vf=2G")
    return stage4_dir


def filt(stage4_dir, dindel_dir):
    raw_indels = getFiles.get_files(stage4_dir, "variantCalls.vcf")
    shell_address = dindel_dir + "/filter_indel.sh"
    SHELL = open(shell_address, 'w')
    SHELL.write("%s\n" % (SGEHEAD))
    for raw_indel in raw_indels:
        indel = dindel_dir + "/" + os.path.basename(raw_indel).split(".")[0]\
                + "_indel.xls"
        SHELL.write("python %s %s %s" % (FILTER, raw_indel, indel))
    SHELL.close()
    manJob.manJobs([shell_address], "-cwd -l vf=1G")


def main(sort_dir, ref_dir, dindel_dir):
    dirCheck.mkdir(dindel_dir)
    tem_dir = dindel_dir + "/tem"
    dirCheck.mkdir(tem_dir)
    stage1_dir = stage1(sort_dir, ref_dir, tem_dir)
    stage2_dir = stage2(stage1_dir, tem_dir)
    stage3_dir = stage3(sort_dir, ref_dir, stage1_dir, stage2_dir, tem_dir)
    stage4_dir = stage4(ref_dir, stage3_dir, tem_dir)
    filt(stage4_dir, dindel_dir)


if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
else:
    print >>sys.stderr, "python %s <sort_dir> <ref_dir> <indel_dir>" %\
    (sys.argv[0])
