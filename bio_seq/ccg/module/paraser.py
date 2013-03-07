#! /usr/bin/python

def paraser(argvs, par_dic, **argv):
    def end_sym(t_str):
        if t_str[0] == '-':
            return False
        for i in xrange(len(t_str)):
            if i > 0:
                if t_str[i-1:i+1] == ' -':
                    return t_str[0:i-1]# Return before next option
        return t_str

    for i in xrange(len(argvs)):
        if argvs[i] in par_dic:
            try:
                par_dic[argvs[i]] = end_sym(' '.join(argvs[i+1:])).split()
            except:
                par_dic[argvs[i]] = True

    if "ARGV_LEN" in par_dic:
        print "[ERROR 1] initail words can't be define!"
    else:
        par_dic["ARGV_LEN"] = len(par_dic) - [par_dic[x] for x in par_dic].count(False)

    return par_dic
