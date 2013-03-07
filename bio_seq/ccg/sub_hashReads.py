#! /usr/bin/python

import sys
import time
import sub_functions

K = 0
F = sub_functions
S = F.ign_direc# To generate a kemr_str whiout forward-relation for hash.

def hash1Read(record, K_graph, K_edges):
    def push2K_graph(kmer, points):
        if kmer not in K_graph:# inital a kmer
            K_graph[kmer] = {'vector':{}, 'libs':{}, 'depth':0}
        kmer_stat = K_graph[kmer]
        for vector in points:
            if vector in kmer_stat['vector']:
                kmer_stat['vector'][vector] += 1
            else:
                kmer_stat['vector'][vector] = 1
        if len(points) == 2:
            K_graph[kmer]['depth'] += 1

    SEQ = record[1]
    length = len(SEQ)
    loop_end = length - K

    if length < K:# no use
        return 1

    for i in xrange(0, loop_end + 1):
        kmer = S(SEQ[i:i+K])

        if i == 0:
            points = (S(SEQ[i+1:i+1+K]),)
        elif i == loop_end:
            points = (S(SEQ[i-1:i-1+K]),)
        else:
            points = (S(SEQ[i-1:i-1+K]), S(SEQ[i+1:i+1+K]))

        push2K_graph(kmer, points)# append kmer and vector(s) to graph

def pe_built(record_1, record_2, insertsize, K_graph, K_edges):
    def lib_push(kmer1, kmer2, lib_name):
        kmer1_stat = K_graph[kmer1]['libs']
        kmer2_stat = K_graph[kmer2]['libs']
        if lib_name not in kmer1_stat:
            kmer1_stat[lib_name] = {}
        if lib_name not in kmer2_stat:
            kmer2_stat[lib_name] = {}

        if kmer2 in kmer1_stat[lib_name]:
            kmer1_stat[lib_name][kmer2] += 1
        else:
            kmer1_stat[lib_name][kmer2] = 1

        if kmer1 in kmer2_stat[lib_name]:
            kmer2_stat[lib_name][kmer1] += 1
        else:
            kmer2_stat[lib_name][kmer1] = 1

    (SEQ1, SEQ2) = record_1[1], record_2[1]
    (len1, len2) = len(SEQ1), len(SEQ2)
    loop_len = min(len1, len2) - K + 1# chose the min read-length to loop
    lib_name = "lib_%d" % insertsize
    #print length

    for i in xrange(loop_len):
        kmer2_endidx = len2 - i

        kmer1 = S(SEQ1[i:i+K])# forward pick from seq_start
        kmer2 = S(SEQ2[kmer2_endidx-K:kmer2_endidx])# forward pick from seq_end

        lib_push(kmer1, kmer2, lib_name)

def remove_low_deepth_edge(K_graph, K_edges, max_multiple = 10):
    def next_seq(kmer_stat):
        best_one = ''
        if 'in' in kmer_stat:
            best_one = chose_best_one(kmer_stat['in'], seq[-1])
        if best_one == '':
            if 'out' in kmer_stat:
                #print kmer_stat
                best_one = chose_best_one(kmer_stat['out'], seq[-1])
        return best_one

    delete_edge = []
    for start_edge in K_edges:
        seq = [chose_start(start_edge, K_graph[start_edge])]
        if seq == ['']:
            continue
        #print "\n%s" % start_edge

        depth_lis = []
        average_depth = 0.0

        while True:
            best_one = next_seq(K_graph[F.ign_direc(seq[-1])])
            if best_one == '':
                #print K_graph[F.ign_direc(seq[-1])]
                delete_edge.append(start_edge)
                break
            depth = node_depth(K_graph[F.ign_direc(best_one)])
            if average_depth == 0:
                average_depth = float(depth)
            else:
                average_depth = sum(depth_lis)/float(len(depth_lis))

            depth_diff_rate = depth / average_depth# can be faster if accept float type error

            if depth_diff_rate > max_multiple:# find error edge
                delete_edge.append(start_edge)
                #F.debug("\n")
                #print K_graph[F.ign_direc(best_one)]
                break

            #F.debug("%d-%.1f," % (depth, depth_diff_rate))
            if depth_diff_rate > 1.5:# normal exit loop
                #F.debug("\n")
                #print K_graph[F.ign_direc(best_one)]
                break

            depth_lis.append(depth)
            seq.append(best_one)

    if len(delete_edge):
        for key in delete_edge:
            del K_edges[key]
            #print key

    return len(delete_edge)

def hashReads(argv):
    hash_time = time.clock()

    global K
    K = argv['kmer']
    K_graph = {}
    K_edges = {'org':{}, 'input':{}, 'filt':{}}
    index = 0

    for fqFile_lis, insertsize in argv['fqFile']:#travel fq libs
        fqFile_len = len(fqFile_lis)

        if fqFile_len == 1:# single end library
            pass
        elif fqFile_len == 2:# pair end library
            (fq_1, fq_2) = F.openFq(fqFile_lis[0]), F.openFq(fqFile_lis[1])

            for record_1 in fq_1:
                index += 1
                record_2 = fq_2.next()

                hash1Read(record_1, K_graph, K_edges)
                hash1Read(record_2, K_graph, K_edges)
                pe_built(record_1, record_2, insertsize, K_graph, K_edges)
        else:
            print >> sys.stderr, "[LIB_ERROR]", fqFile_lis

    depths = [K_graph[kmer]['depth'] for kmer in K_graph]
    min_depth = min(depths)
    avg_depth = sum(depths) / float(len(depths))
    argv['kmer_avg_depth'] = avg_depth

    for kmer in K_graph:
        depth = K_graph[kmer]['depth']
        if depth == 0:
            K_edges['org'][kmer] = 1
        elif depth == min_depth:
            K_edges['input'][kmer] = 1

    hash_time = time.clock() - hash_time
    F.debug("[HASH READS] %d (pe)reads, %d nodes, %d org_edges, %d input_edge, kmer_avg_depth %.2f, %.2fs time_speed\n" % (
        index, len(K_graph), len(K_edges['org']), len(K_edges['input']), avg_depth, hash_time))

    return K_graph, K_edges
