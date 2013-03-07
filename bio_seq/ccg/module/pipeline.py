#! /usr/bin/python

import os
import re

ORG_FUN = {}
ORG_FUN['open'] = open

def debug(msg):
    import sys
    sys.stderr.write(msg + '\n'); sys.exit(0)

def makedir(directory):
    directory = directory.rstrip(os.sep)
    if directory == '':
        return False
    offset = len(directory)
    directories = []
    while offset != -1:
        directories.append(directory[:offset])
        offset = directory[:offset].rfind(os.sep)
    for sub_dir in directories[::-1]:
        if not os.path.exists(sub_dir):
            os.mkdir(sub_dir)
    return True

def open(name, mode='r', **argv):
    if mode == 'w':
        makedir(os.path.dirname(name))# Make sure exist the directory

    if re.compile(r'gz$').match(name):
        import gzip
        return gzip.open(name, mode)
    else:
        return ORG_FUN['open'](name, mode)

def listdir(directory, **argv):
    def strdiff(str_1, str_2):# Single char diff check. if not return False
        if len(str_1) == len(str_2):
            mismatch = []
            for i in xrange(len(str_1)):
                if str_1[i] != str_2[i]:
                    mismatch.append((str_1[i], str_2[i]))
                if len(mismatch)> 1:
                    return False
            if len(mismatch) == 1:
                return mismatch[0]
        return False

    default = {
            'MODE':"NORMAL",
            'SUFFIX':'',
            'DIFF':re.compile(r'^[12]$'),# for "PAIR" module
            }
    for i in default:
        if i not in argv:
            argv[i] = default[i]

    P_suffix = eval("re.compile(r'.*%s$')" % argv['SUFFIX'])
    file_lis = []

    for i in os.listdir(directory):# filt file_lis whatever
        if P_suffix.match(i):
            file_lis.append(i)

    if argv['MODE'] == 'NORMAL':
        pass
    elif argv['MODE'] == 'MANUAL' and argv.has_key('RE'):
        new_lis = []
        P_manual = eval("re.compile(r'%s')" % argv['RE'])
        for i in file_lis:
            if P_manual.match(i):
                new_lis.append(i)
        file_lis = new_lis
    elif argv['MODE'] == 'PAIR':# Find pair fileName, return like [(1.fq, 2.fq), a.fq]
        new_lis = []
        if len(file_lis) >= 2:
            file_lis.sort()

            FIND_PAIR = "Initail"# Iniatal sign variable, except 'False'
            while FIND_PAIR != False:
                if type(FIND_PAIR) == type((0,0)):# Remove the paired record from old list
                    del file_lis[max(FIND_PAIR)]# Fisrt delet max one
                    del file_lis[min(FIND_PAIR)]
                FIND_PAIR = False
                for i in xrange(len(file_lis)):# Traversal finding pair."for 1"
                    if FIND_PAIR != False:# Use the sign variable to break from "for 1"
                        break
                    for j in xrange(len(file_lis)):# "for 2"
                        if j == i:
                            continue
                        diff = strdiff(file_lis[i], file_lis[j])
                        if diff != False:
                            (diff_1, diff_2) = diff
                            if argv['DIFF'].match(diff_1) and argv['DIFF'].match(diff_2):
                                new_lis.append((file_lis[i], file_lis[j]))# Add parired record to new list
                                FIND_PAIR = (i, j)
                                break# break from "for 2"

            file_lis = new_lis + file_lis# Combine the new_list and old_list

    return file_lis

def sameType(obj_input, obj_type):
    if obj_type == type(0):
        try:
            return int(obj_input)
        except:
            debug("[ERROR] %s is not a int num!" % obj_input)
    elif obj_type == type(0.0):
        try:
            return float(obj_input)
        except:
            debug("[ERROR] %s is not a float num!" % obj_input)
    elif obj_type == type(''):
        try:
            return str(obj_input)
        except:
            debug("[ERROR] %s Unkown error, please fix this bug." % obj_input)
    else:
        debug("[ERROR] %s %s Unkown error, please fix this bug." % (
            obj_input, obj_type))

def parser(input_parser, define_parser, **argv):
    def format(parameter, dest):
        sign = dest[0]
        parameter_type = type(define_parser[sign])
        if type(define_parser[sign]) != type([]) and len(dest) > 1:
            debug("[ERROR] '%s' unneed mutiple parameters!" % sign)
        elif type(define_parser[sign]) == type([]):
            parameter_type = type('')
        parameter = sameType(parameter, parameter_type)
        return parameter
    default_argv = {
            'P_sign':re.compile(r'-[\w]$|--[\w]{2,}$'),
            'HELP_CMD':['-h', '--help'],
            'USAGE':'Usage:\n   Please see the program...',
            }
    for i in default_argv:
        if i not in argv:
            argv[i] = default_argv[i]
    for i in define_parser:
        if i in argv['HELP_CMD']:
            debug("[PROAM ERROR] '%s' can't be used as parameter sign."
                    " program default: HELP_CMD=['-h', '--help']" %
            ' '.join(argv['HELP_CMD']))

    pars = []
    if len(input_parser) == 0:
        debug(argv['USAGE'])
    if input_parser[0][0] != '-':
        if 'DEFAULT' in argv:
            pars.append([argv['DEFAULT']])
        else:
            debug(argv['USAGE'])
    for i in input_parser:
        if i == '':
            debug("[ERROR] illegal input !")
        elif i[0] == '-':
            if argv['P_sign'].match(i) == None:
                debug("[ERROR] '%s' :illegal parameter! exp:-w --word" % i)
            elif i in argv['HELP_CMD']:
                debug(argv['USAGE'])
            elif i not in define_parser:
                debug("[ERROR] '%s' :illegal input!" % i)
            else:
                pars.append([i])
        else:
            pars[-1].append(format(i, pars[-1]))

    for i in pars:# Add the input to default
        if type(define_parser[i[0]]) == type([]):
            define_parser[i[0]] = i[1:]
        elif type(define_parser[i[0]]) == type(True):
            define_parser[i[0]] = True
        else:
            define_parser[i[0]] = i[1]

    for i in define_parser:# Check
        if define_parser[i] == '' or define_parser[i] == []:
            debug("[ERROR] '%s' must have parameter(s)." % i)

    return define_parser
