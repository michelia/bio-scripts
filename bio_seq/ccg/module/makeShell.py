#! /usr/bin/python
import os

HEADER = "#$ -S /bin/sh\n#! /bin/bash\n"
HEADER_LEN = len(HEADER)

def formatCon(shell_con):
    try:
        if shell_con[0:HEADER_LEN] != HEADER:
            shell_con = HEADER + shell_con
    except:
        shell_con = HEADER + shell_con

    return shell_con

def makeShell(shell_dir, shell_name, shell_con):#output shell_file, return address.

    shell_con = formatCon(shell_con)
    sep = os.sep
    shell_dir = shell_dir.rstrip(os.sep)
    shell_address = "%s%s%s" % (shell_dir, sep, shell_name)

    OUT = open(shell_address, 'w')
    OUT.write("%s\n" % shell_con)
    OUT.close()

    return shell_address

def checkShell(shell_address):
    formated = []
    for add in shell_address:
        SH_FILE = open(add, 'r')
        content = ''.join(SH_FILE.readlines())
        SH_FILE.close()

        contentFmt = formatCon(content)
        if len(content) != len(contentFmt):
            OUT = open(add, 'w')
            OUT.write(contentFmt)
            OUT.close()
            formated.append(add)

    if len(formated) > 0:
        return formated
    else:
        return True
