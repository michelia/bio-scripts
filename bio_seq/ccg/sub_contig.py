#! /usr/bin/python

import sys
import time
import sub_functions

F = sub_functions
S = F.ign_direc

def chose_start(K_graph, this_edge, visited, argv, mode = 'seq'):
    #print this_edge
    for visite_kmer in visited:
        if visite_kmer in K_graph[this_edge]['vector']:
            del K_graph[this_edge]['vector'][visite_kmer]

    edge_next = K_graph[this_edge]['vector']
    (start_1, start_2) = this_edge.split('|')
    start_lis = []
    next_lis = []

    for vector in edge_next:
        (fmt1, fmt2) = vector.split('|')

        if fmt1.find(start_1[1:]) == 0:
            start_lis.append(start_1)
            next_lis.append(fmt1)
        elif fmt2.find(start_1[1:]) == 0:
            start_lis.append(start_1)
            next_lis.append(fmt2)

        if fmt1.find(start_2[1:]) == 0:
            start_lis.append(start_2)
            next_lis.append(fmt1)
        elif fmt2.find(start_2[1:]) == 0:
            start_lis.append(start_2)
            next_lis.append(fmt2)

    start_lis_len = len(start_lis)
    next_lis_len = len(next_lis)

    if start_lis_len == 0 and next_lis_len == 0:
        return []
    elif start_lis_len == 1 and next_lis_len == 1:
        return start_lis
    elif start_lis_len == 2 and next_lis_len == 2:
        if mode == 'seq':
            del K_graph[this_edge]['vector'][S(next_lis[1])]
        res = [start_lis[0], next_lis[1]]
        return res
    else:
        print >> sys.stderr, start_lis, next_lis
        print >> sys.stderr, "[DEBUG] chose_start error.", K_graph[this_edge]; sys.exit()
        return False

def seqlis2str(seq_lis):
    res = [seq_lis[0]]
    for i in xrange(1, len(seq_lis)):
        res.append(seq_lis[i][-1])
    res = ''.join(res)
    return res

def seqRstrip(seq_lis):
    del_seq = seq_lis[-1]
    while seq_lis[-2] == del_seq:
        del seq_lis[-1]
    return seq_lis

def checkIsCircle(seq, artif_seq_lis):
    for i in xrange(0, len(artif_seq_lis) - 1):
        if seq[-i] != artif_seq_lis[-i]:
            return False
    print >> sys.stderr, "[CIRCLE] have circle"
    return True

def nodes2contig(K_graph, edges_lis, visited, artif_seq_lis, argv):
    global K_index
    true_rate = argv['true_rate']
    if len(artif_seq_lis) > 1:
        mode = 'extend'
        base_seq = seqlis2str(artif_seq_lis)
    else:
        mode = 'contig'

    if mode == 'extend':# extend mode
        seq = artif_seq_lis + []
    elif mode == 'contig':#contig mode
        this_edge = edges_lis[0]# pick one to extend
        del edges_lis[0]# used edge

        if this_edge in visited:# visited edge
            return False
        else:
            visited[this_edge] = 1# push to visited_dic because used for begining

        start_seq = chose_start(K_graph, this_edge, visited, argv)# will delete used kmer in start_kemr_vectors
        start_seq_len = len(start_seq)

        if start_seq_len == 1:
            seq = [start_seq[0]]
        elif start_seq_len == 0:
            print >> sys.stderr, "[ERROR] start edge error", this_edge, edges_lis
            return False
        else:
            edges_lis.extend([S(x) for x in start_seq])# split a two-direction edge
            #print "!", edges_lis
            return False

    #F.debug("[DEBUG] %s donging %s..." % (seq[-1], mode))

    while 1:
        next_seq = F.chose_best_one(K_graph, seq, argv)
        next_seq_len = len(next_seq)

        if next_seq_len == 1:# unique result return
            seq.append(next_seq[0])
            if seq[-1] == seq[-2]:# may duplicate
                min_insert_size = argv['insert_size'][0][0]
                if len(seq) > min_insert_size:
                    if seq[-1] == seq[min_insert_size]:
                        seq = seqRstrip(seq)
                        visited[S(seq[-1])] = 1
                        break

            #jump form circle extend loop
            if mode == 'extend':
                if len(seq) > len(artif_seq_lis):
                    if checkIsCircle(seq, artif_seq_lis):
                        break
            #print len(seqlis2str(seq))

        elif next_seq_len > 1:
            if mode == 'contig':
                edges_lis.extend([S(x) for x in next_seq])
                visited[S(seq[-1])] = 1# push to the visited because use for terminal
                break
            elif mode == 'extend':
                break
        else:# unusual event
            #print seq[-1], K_graph[S(seq[-1])]; sys.exit()
            visited[S(seq[-1])] = 1
            break

    con = seqlis2str(seq)
    K_index += 1

    return (">contig%d_k%d_len%d" % (K_index, K, len(con)), con)# return name, seq_str

def removeLowEdges(K_graph, K_edges, argv):
    def isLowQualEdge(edge, true_rate):
        start = chose_start(K_graph, edge, {}, argv, 'check')
        start_len = len(start)
        if start_len != 1:
            return False

        seq = [start[0]]
        depths = []

        for i in xrange(argv['read_len']):
            this_edge = S(seq[-1])
            vector_dic = K_graph[this_edge]['vector']
            next_lis = F.extend_tbl(seq[-1], vector_dic)
            next_lis_len = len(next_lis)

            if next_lis_len == 1:
                seq.append(next_lis[0][0])
                depths.append(K_graph[this_edge]['depth'])
            elif next_lis_len == 0:
                return True
            else:
                if next_lis[0][1] >= true_rate:
                    seq.append(next_lis[0][0])
                    depths.append(K_graph[this_edge]['vector'][S(seq[-1])])

            depths_len = len(depths) - 1
            depths_sum =float( sum(depths[0:-1]))
            if depths_len > 1:
                avg_depth = depths_sum / depths_len
                diff_rate = 1 - (avg_depth / depths[-1])
                if diff_rate >= true_rate:
                    #print diff_rate, depths, edge# debug: false edge trend
                    return True

        #print depths# debug: right edge trend
        return False

    true_rate = argv['true_rate']
    true_rate = 0.9# debug
    need_del_org = []
    need_del_input = []

    for edge in K_edges['org']:
        if isLowQualEdge(edge, true_rate):
            need_del_org.append(edge)
    for edge in need_del_org:
        del K_edges['org'][edge]

    for edge in K_edges['input']:
        if isLowQualEdge(edge, true_rate):
            need_del_input.append(edge)
    for edge in need_del_input:
        del K_edges['input']

    #print K_edges; print len(need_del_org), need_del_org# debug
    return len(need_del_org) + len(need_del_input)# return deleted edge(s) num

def graph2contig(K_graph, K_edges, argv):

    argv['insert_size'] = sorted(argv['insert_size'], key = lambda x:x[0])
    global K, K_index
    K = argv['kmer']
    K_index = 0

    remove_time = time.clock()
    remove_num = removeLowEdges(K_graph, K_edges, argv)
    remove_time = time.clock() - remove_time
    F.debug("[EDGES] remove %d low quality edge(s), %d remain, %.2f time_speed\n" % (remove_num,
        len(K_edges['org']) + len(K_edges['input']), remove_time))

    kmer_edges = [x for x in K_edges['org']] + [x for x in K_edges['input']]# combine edges
    K_edges['visited'] = {}

    contigs = []
    contig_time = time.clock()

    while len(kmer_edges):
        contig = nodes2contig(K_graph, kmer_edges, K_edges['visited'], [], argv)
        if contig != False:
            #print len(contig[1])# debug
            contigs.append(contig)

    contig_time = time.clock() - contig_time
    contig_len = sum([len(x[1]) for x in contigs])
    F.debug("[CONTIG] %d contigs, %d total_len, %.2f time_speed\n" % (len(contigs),
        contig_len, contig_time))

    return contigs

def extendFasta(K_graph, K_edges, argv):
    def seq2seq_lis(seq_str):
        seq_str_len = len(seq_str)
        res = []
        for i in xrange(seq_str_len):
            if i + K > seq_str_len:
                break
            res.append(seq_str[i:i+K])
        return res
    def combine(f_seq_str, r_seq_str, ref_len):
        r_seq_str = F.reverse(r_seq_str)
        if r_seq_str[-ref_len:] == f_seq_str[0:ref_len]:
            return r_seq_str + f_seq_str[ref_len:]
        else:
            print "[ERROR] extend error" >> sys.stderr
            sys.exit()

    argv['insert_size'] = sorted(argv['insert_size'], key = lambda x:x[0])
    global K, K_index
    K = argv['kmer']
    K_index = 0

    #remove_time = time.clock()
    #remove_num = removeLowEdges(K_graph, K_edges, argv)
    #remove_time = time.clock() - remove_time
    #F.debug("[EDGES] remove %d low quality edge(s), %d remain, %.2f time_speed\n" % (remove_num,
    #    len(K_edges['org']) + len(K_edges['input']), remove_time))

    seqs = argv['refSeq']
    contigs = []
    contig_time = time.clock()

    for name, seq in seqs:
        seq_len = len(seq)
        if seq_len <= K:
            contigs.append((name, seq))
            continue
        seq_r = F.reverse(seq)
        seq_ext_lis = nodes2contig(K_graph, [], {}, seq2seq_lis(seq), argv)
        seq_r_ext_lis = nodes2contig(K_graph, [], {}, seq2seq_lis(seq_r), argv)
        seq_extended = combine(seq_ext_lis[1], seq_r_ext_lis[1], seq_len)
        contigs.append(("%s old_len %d new_len %d" % (name, len(seq), len(seq_extended)), seq_extended))

    contig_time = time.clock() - contig_time
    contig_len = sum([len(x[1]) for x in contigs])
    F.debug("[CONTIG] %d contigs, %d total_len, %.2f time_speed\n" % (len(contigs),
        contig_len, contig_time))

    return contigs
