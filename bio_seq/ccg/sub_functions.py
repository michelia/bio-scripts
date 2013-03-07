#! /usr/bin/python

import sys

def debug(msg, exit = False):
    sys.stderr.write(msg)
    if exit != False:
        sys.exit(exit)

def reverse(seq):
    dic = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    res = []
    for i in seq:
        res.append(dic[i])
    return ''.join(res[::-1])

def ign_direc(kmer):
    return '|'.join(sorted((kmer,reverse(kmer))))

def openFq(fqFile):
    IN = open(fqFile)
    next = IN.next
    while 1:
        lines = (next().rstrip(), next().rstrip(), next().rstrip(), next().rstrip())
        yield lines
    IN.close()

def extend_tbl(seq, vector_dic):
    res = []
    match_t = 0.0
    for vector in vector_dic:
        (fmt1, fmt2) = vector.split('|')
        if fmt1.find(seq[1:]) == 0:
            res.append((fmt1, vector_dic[vector]))
            match_t += vector_dic[vector]
        elif fmt2.find(seq[1:]) == 0:
            res.append((fmt2, vector_dic[vector]))
            match_t += vector_dic[vector]

    if match_t == 0:
        #print seq, vector_dic; sys.exit()
        return res
    res = [(x[0], x[1]/match_t) for x in res]# compute prob rate
    res = sorted(res, key = lambda x:x[1], reverse=True)# sort

    return res

def libRelation_num(K_graph, seq, best_one_lis, libs_lis, argv):
    #print extend_tbl_lis
    seq_len = len(seq)
    rel_vector_dic = {}
    filter_num = 0
    index = 0
    libRelation_tbl_lis = []

    for next_seq in best_one_lis:
        rel_vector_dic[ign_direc(next_seq)] = 0

    for lib_size, sd in libs_lis:
        read_len = argv['read_len']
        if False:# whether use this lib for relation fix
            continue

        lib_name = "lib_%d" % lib_size
        back_point = seq_len - lib_size + argv['read_len']
        back_fr = back_point - int(1.5 * sd)
        back_to = back_point + int(1.5 * sd)
        back_seq_dic = {}

        for i in xrange(back_fr, back_to + 1):# base on back back-kmers, scaf issue-kmer
            seq_i_kmer = ign_direc(seq[i])
            back_seq_dic[seq_i_kmer] = 1
            for kmer in K_graph[seq_i_kmer]['libs'][lib_name]:
                index += 1
                if kmer in rel_vector_dic:
                    filter_num += 1
                    rel_vector_dic[kmer] += 1

    next_seq_num = []
    for next_seq in best_one_lis:
        next_seq_num.append([next_seq, rel_vector_dic[ign_direc(next_seq)]])

    return next_seq_num

def libRelation_tbl(K_graph, seq, extend_tbl_lis, libs_lis, argv):
    #print extend_tbl_lis
    seq_len = len(seq)
    rel_vector_dic = {}
    filter_num = 0
    index = 0
    libRelation_tbl_lis = []

    for next_seq, prob in extend_tbl_lis:
        rel_vector_dic[ign_direc(next_seq)] = 0.0

    for lib_size, sd in libs_lis:
        read_len = argv['read_len']
        if False:# whether use this lib for relation fix
            continue

        lib_name = "lib_%d" % lib_size
        back_point = seq_len - lib_size + argv['read_len']
        back_fr = back_point - int(2 * sd)
        back_to = back_point + int(2 * sd)
        back_seq_dic = {}

        for i in xrange(back_fr, back_to + 1):# base on back back-kmers, scaf issue-kmer
            seq_i_kmer = ign_direc(seq[i])
            back_seq_dic[seq_i_kmer] = 1

            if lib_name not in K_graph[seq_i_kmer]['libs']:
                continue

            for kmer in K_graph[seq_i_kmer]['libs'][lib_name]:
                index += 1
                if kmer in rel_vector_dic:
                    filter_num += 1
                    rel_vector_dic[kmer] += 1

        #for next_seq, prob in extend_tbl_lis:# base on issue-kmer, scaf back-kmers
        #    next_seq_kmer = ign_direc(next_seq)
        #    for vector in K_graph[next_seq_kmer]['libs'][lib_name]:
        #        if vector in back_seq_dic:
        #            rel_vector_dic[next_seq_kmer] += 1

    #print index, rel_vector_dic, seq_len
    if filter_num == 0:
        return extend_tbl_lis

    total_prob = 0.0
    for next_seq, prob in extend_tbl_lis:
        kmer = ign_direc(next_seq)
        prob = prob * (rel_vector_dic[kmer]/filter_num)
        total_prob += prob
        if prob == 0:
            continue
        libRelation_tbl_lis.append((next_seq, prob))

    libRelation_tbl_lis = [(x[0], x[1]/total_prob) for x in libRelation_tbl_lis]
    return libRelation_tbl_lis

def combineNextLis(combine_next, combine_next_base_pro):
    ################ [chosed-base, [[base-seq, prob, relaton_kmer_num]]]
    #print len(combine_next)
    total_match = 0.0

    for j in xrange(len(combine_next)):# static total_match
        (next_seq, seq_pr_lis) = combine_next[j]
        for k in xrange(len(seq_pr_lis)):# every batch
            (seq_pr, prob, kmer_nums) = seq_pr_lis[k]
            #total_match += kmer_nums
            total_match += kmer_nums * combine_next_base_pro[j]

    for j in xrange(len(combine_next)):# compute prob
        (next_seq, seq_pr_lis) = combine_next[j]
        for k in xrange(len(seq_pr_lis)):
            (seq_pr, prob, kmer_nums) = seq_pr_lis[k]

            if total_match == 0:
                prob = 0
            else:
                #prob = kmer_nums / total_match
                prob = (kmer_nums * combine_next_base_pro[j])/ total_match

            #print next_seq, seq_pr[-1], prob
            combine_next[j][1][k] = [seq_pr, prob, kmer_nums]
            #print len(combine_next[0][1][0])

    for j in xrange(len(combine_next)):# sort batch by prob
        (next_seq, seq_pr_lis) = combine_next[j]
        tmp = seq_pr_lis
        tmp = sorted(tmp, key = lambda x:x[1], reverse = True)
        combine_next[j][1] = sorted(seq_pr_lis, key = lambda x:x[1], reverse = True)

    combine_next = sorted(combine_next, key = lambda x:x[1][0][1], reverse = True)# sort stcik by prob

    #print total_match
    #print combine_next[0][1][0][1]
    return combine_next

def chose_best_one(K_graph, seq, argv):
    #print argv
    #libs_lis = sorted(argv['insert_size'], key = lambda x:x[0])
    libs_lis = argv['insert_size']
    true_rate = argv['true_rate']
    trans_rate = argv['trans_rate']
    vector_dic = K_graph[ign_direc(seq[-1])]['vector']

    #print seq[-1]
    seq_len = len(seq)
    next_lis =  extend_tbl(seq[-1], vector_dic)
    next_lis_len = len(next_lis)

    if next_lis_len == 1:
        return [next_lis[0][0]]
    elif next_lis_len == 0 :
        return []
    else:
        if next_lis[0][1] >= true_rate:# prob accept
            return [next_lis[0][0]]

        if len(libs_lis) == 0:# no libs
            return [x[0] for x in next_lis]

        (min_ins, min_ins_sd) = libs_lis[0]

        if len(seq) < min_ins - argv['read_len'] + min_ins_sd * 1.5:# too short to analy
            return [x[0] for x in next_lis]

        ##################### libRelation correted ########################
        lib_next_lis = libRelation_tbl(K_graph, seq, next_lis, libs_lis, argv)

        if len(lib_next_lis) == 1:
            return [lib_next_lis[0][0]]
        elif len(lib_next_lis) == 0:
            return []
        else:# > 1
            if lib_next_lis[0][1] >= true_rate:
                return [lib_next_lis[0][0]]

            #print next_lis, lib_next_lis# debug

            ################ [chosed-base, [[base-seq, prob, relaton_kmer_num]]]
            combine_next = [ [x[0], [[seq + [x[0]], x[1], 0]]] for x in lib_next_lis ]

            max_try_times = libs_lis[-1][0]# max lib_width
            for i in xrange(max_try_times):

                try_combine_next = []
                try_combine_next_base_pro = []
                for j in xrange(len(combine_next)):# every prob-base
                    (next_seq, seq_pr_lis) = combine_next[j]
                    seq_pr_lis_new = []

                    for k in xrange(len(seq_pr_lis)):# every batch
                        (seq_pr, prob, kmer_nums) = seq_pr_lis[k]

                        best_one_lis = chose_best_one(K_graph, seq_pr, argv)
                        best_one_lis_len = len(best_one_lis)
                        best_one_lis_rls = libRelation_num(K_graph, seq_pr, best_one_lis, libs_lis, argv)

                        for best_one, rel_num in best_one_lis_rls:
                            seq_pr_lis_new.append([seq_pr + [best_one], 0, rel_num])

                    if len(seq_pr_lis_new):
                        try_combine_next.append([next_seq, seq_pr_lis_new])# if seq_pr_lis_new == [], this maybe extend is termianal
                        try_combine_next_base_pro.append(prob)

                if len(try_combine_next) == 0:# Unknow situation, return None
                    return []

                combine_next = combineNextLis(try_combine_next, try_combine_next_base_pro)# statical sort

                max_prob = combine_next[0][1][0][1]
                max_prob_next_seq = combine_next[0][0]

                if max_prob == 0:
                    return []

                if max_prob >= true_rate:
                    #print max_prob, max_prob_next_seq
                    #sys.exit()
                    return [max_prob_next_seq]# return prob-base
                #print max_prob, max_prob_next_seq

            return []# can not fix with lib relation

