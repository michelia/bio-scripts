#! /usr/bin/python

import sys
import re
import sys
import math

def qualitySta(quality):
    qualitys = []
    for i in range(len(quality)):
        per_quality = int(ord(quality[i])) - 33
        qualitys.append(per_quality)
    return qualitys


def avgQuality(qualitys):
    qua_num = len(qualitys)
    avg_qua = 0
    if qua_num > 0:
        total_qua = 0
        for qua in qualitys:
            total_qua += qua
        avg_qua = total_qua / qua_num
    return avg_qua


def consensusMatchAnalysis(match):
    matchs = []
    largest_indel = 0
    while True:
        add_len = 1
        alt = match[0]
        pre_alt = ''
        if len(match) > 1:
            pre_alt = match[1]
        if alt == '+' or alt == '-':  #indel
            if ord(match[1]) >= 48 and ord(match[1]) <= 57:
                if int(match[1]) > largest_indel:
                    largest_indel = int(match[1])
                add_len = int(match[1]) + 2
                if ord(match[2]) >= 48 and ord(match[2]) <= 57:
                    add_len = int(match[1:3]) + 3
                    if int(match[1:3]) > largest_indel:
                        largest_indel = int(match[1:3])
        elif alt == '^':  #begin
            add_len = 3
        elif pre_alt == '$':  #end
            add_len = 2
        alt = match[:add_len]
        match = match[add_len:]

        if not re.search('^[+-]', alt):
            matchs.append(alt)
        if match == '':
            break

    return (matchs, largest_indel)


def consensusAnalysis(line, read_len):

    (chrom, coordinate, reference, depth, match, baq_quality, map_quality,\
            position) = line.split()

    (matchs, indel_size) = consensusMatchAnalysis(match)
    map_qualitys = qualitySta(map_quality)
    baq_qualitys = qualitySta(baq_quality)
    positions = position.split(',')

    if len(matchs) != len(map_qualitys) or len(matchs) != len(positions) or\
            len(matchs) != len(baq_qualitys):  #check above analysis
        print >>sys.stderr, "statistics error at coordinate: ", coordinate
        sys.exit(1)

    match_dic = {'A' : [0, 0, 0, [], []],
            'T' : [0, 0, 0, [], []],
            'G' : [0, 0, 0, [], []],
            'C' : [0, 0, 0, [], []]
            }  #total, forward, reverse, baqs, maqs

    for i in range(len(matchs)):
        per_match = matchs[i]
        per_maq_qua = map_qualitys[i]
        per_baq_qua = baq_qualitys[i]
        per_position = int(positions[i])

        maq_lower = 0
        if per_position < 5 or per_position + 4 > read_len:
            maq_lower = -10

        match_base = ''
        if len(per_match) == 1:
            match_base = per_match

        elif re.search('^\^', per_match):
            match_base = per_match[-1]
            maq_lower = -15

        elif re.search('\$$', per_match):
            match_base = per_match[0]
            maq_lower = -15

        if match_base == '':
            print >>sys.stderr, match
            sys.exit(1)
        else:
            if match_base == '.':
                match_base = reference
            elif match_base == ',':
                match_base = reference.lower()
            elif match_base == 'N' or match_base == 'n':
                continue
            elif match_base == '>' or match_base == '<':
                continue
            elif match_base != '*':
                match_base = match_base
            else:
                continue

            per_maq_qua += maq_lower

            if per_maq_qua < 0:
                per_maq_qua = 0

            if match_base in match_dic:  #forward
                match_dic[match_base][1] += 1
            else:  #reverse
                match_dic[match_base.upper()][2] += 1

            match_dic[match_base.upper()][0] += 1
            match_dic[match_base.upper()][3].append(per_baq_qua)
            match_dic[match_base.upper()][4].append(per_maq_qua)

    return (match_dic, indel_size)


def baseAnalysisX(chrom, coordinate, ref, match_dic, min_baq_avg, min_maq_avg):
    new_matchs = sorted(match_dic.items(), key = lambda d:d[1][0])
    base1 = new_matchs[-1][0]
    base1_depth = new_matchs[-1][1][0]
    base1_depth_forward = new_matchs[-1][1][1]
    base1_depth_reverse = new_matchs[-1][1][2]
    base1_baq_avg = avgQuality(new_matchs[-1][1][3])
    base1_maq_avg = avgQuality(new_matchs[-1][1][4])

    base2 = new_matchs[-2][0]
    base2_depth = new_matchs[-2][1][0]
    base2_depth_forward = new_matchs[-2][1][1]
    base2_depth_reverse = new_matchs[-2][1][2]
    base2_baq_avg = avgQuality(new_matchs[-2][1][3])
    base2_maq_avg = avgQuality(new_matchs[-2][1][4])

    low_quality_info = "LOWQUA\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\
\t%d\t%d\t%d" % (chrom, coordinate, ref, base1, base1_depth,\
            base1_depth_forward, base1_depth_reverse, base1_baq_avg,\
            base1_maq_avg, base2, base2_depth, base2_depth_forward,\
            base2_depth_reverse, base2_baq_avg, base2_maq_avg)

    is_mismatch = 0
    if (base1 != ref and base1_depth != 0) or\
            (base2 != ref and base2_depth != 0):
        is_mismatch = 1

    is_low_quality = 0
    if base1_depth == 0:
        base1 = 'N'
        base1_depth = 0
        base1_depth_forward = 0
        base1_depth_reverse = 0
    if base1_baq_avg < min_baq_avg or base1_maq_avg < min_maq_avg:
        base1 = 'N'
        base1_depth = 0
        base1_depth_forward = 0
        base1_depth_reverse = 0
        is_low_quality = 1

    if base2_depth == 0 or base1_depth == 0:
        base2 = 'N'
        base2_depth = 0
        base2_depth_forward = 0
        base2_depth_reverse = 0
    if base2_baq_avg < min_baq_avg or base2_maq_avg < min_maq_avg:
        base2 = 'N'
        base2_depth = 0
        base2_depth_forward = 0
        base2_depth_reverse = 0
        is_low_quality = 1

    if is_mismatch == 1 and is_low_quality == 1:
        print >>sys.stderr, low_quality_info
    return (base1, base1_depth, base1_depth_forward, base1_depth_reverse,\
            base1_baq_avg, base1_maq_avg, base2, base2_depth,\
            base2_depth_forward, base2_depth_reverse, base2_baq_avg,\
            base2_maq_avg)


def gcPercentage(sequences, read_len):
    gc_up = 0
    for i in range(read_len - 10, read_len):
        if sequences[i] == 'G' or sequences[i] == 'C':
            gc_up += 1

    gc_up_rate = gc_up / 10.0

    gc_down = 0
    for j in range(read_len + 1, read_len + 11):
        if sequences[j] == 'G' or sequences[j] == 'C':
            gc_down += 1
    gc_down_rate = gc_down / 10.0

    if abs(gc_up_rate - 0.5) >= abs(gc_down_rate - 0.5):
        return gc_up_rate
    else:
        return gc_down_rate


def snpEdge(sequences, base1, base2, read_len):
    base1_up = 0
    base1_down = 0
    base2_up = 0
    base2_down = 0
    for i in range(read_len - 1, -1, -1):
        if sequences[i] == base1:
            base1_up += 1
        else:
            break
    for i in range(read_len - 1, -1, -1):
        if sequences[i] == base2:
            base2_up += 1
        else:
            break

    for j in range(read_len + 1, len(sequences)):
        if sequences[j] == base1:
            base1_down += 1
        else:
            break
    for j in range(read_len + 1, len(sequences)):
        if sequences[j] == base2:
            base2_down += 1
        else:
            break
    return (base1_up, base1_down, base2_up, base2_down)


def modDepth(depths, read_length):
    depth_dic = {}
    for i in range(read_length - 10, read_length + 11):
        depth = depths[i]
        if depth in depth_dic:
            depth_dic[depth] += 1
        else:
            depth_dic[depth] = 1
    mod_depth = sorted(depth_dic.items(), key = lambda d:d[1])[-1][0]
    return mod_depth


def mismatchDensity(mismatch_distribution, coordinates, window, read_len):
    #analysis read length position
    #window = read_len * 2 - 1
    mismatch_sta = 0
    shortest_distance = 1000  #means very large
    for i in range(window):
        if mismatch_distribution[i] > 0:  #> 1 maybe better
            mismatch_sta += 1
            difference = abs(mismatch_distribution[i] - mismatch_distribution\
                    [read_len]) / float(mismatch_distribution[i] +\
                    mismatch_distribution[read_len])
            if difference < 0.6:
                distance = abs(coordinates[i] - coordinates[read_len]) - 1
                if distance > -1 and distance < shortest_distance:
                    shortest_distance = distance
    return (mismatch_sta, shortest_distance)


def baseSequencingAnalysis(base1, depth1, depth1_forward, depth1_reverse,\
        avg_baq1, avg_maq1, forward_duplication1, reverse_duplication1,\
        base2, depth2, depth2_forward, depth2_reverse, avg_baq2, avg_maq2,\
        forward_duplication2, reverse_duplication2, error_info):

        error_rate = 0.05  #allert error rate is 0.05 ~ Q13
        depth = depth1 + depth2

        is_sequencing_error1 = 0
        duplication1 = forward_duplication1
        if reverse_duplication1 > duplication1:
            duplication1 = reverse_duplication1
        error_depth1 = depth * error_rate * (duplication1 ** duplication1)
        if error_depth1 >= (depth1 - 1):  #sequencing error 1
            is_sequencing_error1 = 1
        if depth1_forward * depth1_reverse == 0:  #sequencing error 2
            if duplication1 > 0 and avg_baq1 < 26:
                is_sequencing_error1 = 1
        if is_sequencing_error1 == 1:
            base1 = 'N'
            depth1 = 0
            depth1_forward = 0
            depth1_reverse = 0

        is_sequencing_error2 = 0
        duplication2 = forward_duplication2
        if reverse_duplication2 > duplication2:
            duplication2 = reverse_duplication2
        error_depth2 = depth * error_rate * (duplication2 ** duplication2)
        if error_depth2 >= (depth2 - 1):  #sequencing error 1
            is_sequencing_error2 = 1
        if depth2_forward * depth2_reverse == 0:  #sequencing error 2
            if duplication2 > 0 and avg_baq2 < 26:
                is_sequencing_error2 = 1
        if is_sequencing_error2 == 1:
            base2 = 'N'
            depth2 = 0
            depth2_forward = 0
            depth2_reverse = 0

        if is_sequencing_error1 == 1 or is_sequencing_error2 == 1:
            print >>sys.stderr, 'SEQUENCINGERROR', error_info
        return (base1, depth1, depth1_forward, depth1_reverse, base2, depth2,\
                depth2_forward, depth2_reverse)


def baseChange(base1, base1_depth, base1_depth_forward, base1_depth_reverse,\
        base2, base2_depth, base2_depth_forward, base2_depth_reverse,\
        reference):
    is_change = 0
    if base1_depth > 0 and base2_depth > 0:
        if base2 == reference:
            is_change = 1
        elif base1 != reference and base2 != reference and\
                (base2_depth > base1_depth):
            is_change = 1
    elif base1_depth == 0 and base2_depth > 0:
        if base2 != reference:
            base1 = reference
        else:
            print >>sys.stderr, 'base change error 1', base1, base1_depth,\
                    base1_depth_forward, base1_depth_reverse, base2,\
                    base2_depth, base2_depth_forward, base2_depth_reverse,\
                    reference
    elif base1_depth > 0 and base2_depth == 0:
        if base1 != reference:
            base2 = reference
            is_change = 1
        else:
            print >>sys.stderr, 'base change error 2', base1, base1_depth,\
                    base1_depth_forward, base1_depth_reverse, base2,\
                    base2_depth, base2_depth_forward, base2_depth_reverse,\
                    reference

    if is_change == 1:
        tem_base1 = base1
        tem_base1_depth = base1_depth
        tem_base1_depth_forward = base1_depth_forward
        tem_base1_depth_reverse = base1_depth_reverse

        base1 = base2
        base1_depth = base2_depth
        base1_depth_forward = base2_depth_forward
        base1_depth_reverse = base2_depth_reverse

        base2 = tem_base1
        base2_depth = tem_base1_depth
        base2_depth_forward = tem_base1_depth_forward
        base2_depth_reverse = tem_base1_depth_reverse

    return (base1, base1_depth, base1_depth_forward, base1_depth_reverse,\
            base2, base2_depth, base2_depth_forward, base2_depth_reverse)


def baseAnalysisY(pileup, min_baq, min_maq, min_rate, read_len, max_depth):
    snp_for_z = []  #baseAnalysisY result to next step
    matchs = []  #do not analysis low quality match reads now
    coordinates = []
    sequences = []
    mismatch_distribution = []
    depths = []
    window = read_len * 2 - 1
    local_indel_dic = {}

    for line in open(pileup):
        line = line.rstrip()
        ref = line.split()[2]
        chrom = line.split()[0]
        coordinate = int(line.split()[1])

        coordinates.append(coordinate)
        sequences.append(ref)

        (match_dic, indel_size) = consensusAnalysis(line, read_len)

        if indel_size > 0: #local_indel
            for i in range(coordinate, coordinate + indel_size):
                local_indel_dic[i] = 0

        (base1, base1_depth, base1_depth_forward, base1_depth_reverse,\
                base1_baq_avg, base1_maq_avg, base2, base2_depth,\
                base2_depth_forward, base2_depth_reverse, base2_baq_avg,\
                base2_maq_avg) = baseAnalysisX(chrom, coordinate, ref,\
                match_dic, min_baq, min_maq)

        matchs.append((base1, base1_depth, base1_depth_forward,\
                base1_depth_reverse, base1_baq_avg, base1_maq_avg,\
                base2, base2_depth, base2_depth_forward, base2_depth_reverse,\
                base2_baq_avg, base2_maq_avg))

        mismatch = 0
        if base1 != 'N' and base1 != ref:
            mismatch = base1_depth
        elif base2 != 'N' and base2 != ref:
            mismatch = base2_depth
        mismatch_distribution.append(mismatch)

        depths.append(base1_depth + base2_depth)

        if len(coordinates) == window:  #analysis at middle (read len position)

            if mismatch_distribution[read_len] > 0: #find a mismatch
                m_coordinate = coordinates[read_len]
                m_reference = sequences[read_len]
                m_depth = depths[read_len]

                (mismatch_density, shorest_mismatch_distance)= mismatchDensity\
                        (mismatch_distribution, coordinates, window, read_len)

                gc_rate = gcPercentage(sequences, read_len) #find highest gc rate

                (m_base1, m_base1_depth, m_base1_depth_forward,\
                        m_base1_depth_reverse, m_base1_baq_avg,\
                        m_base1_maq_avg, m_base2, m_base2_depth,\
                        m_base2_depth_forward, m_base2_depth_reverse,\
                        m_base2_baq_avg, m_base2_maq_avg) = matchs[read_len]

                (base1_up, base1_down, base2_up, base2_down) =\
                        snpEdge(sequences, m_base1, m_base2, read_len)

                mod_depth = modDepth(depths, read_len)

                #all information to show error snp
                error_snp_info ="%s\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t\
%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%d\t%d" % (chrom, m_coordinate,\
m_reference, m_depth, m_base1, m_base1_depth, m_base1_depth_forward,\
m_base1_depth_reverse, m_base1_baq_avg, m_base1_maq_avg, m_base2,\
m_base2_depth, m_base2_depth_forward, m_base2_depth_reverse, m_base2_baq_avg,\
m_base2_maq_avg, base1_up, base1_down, base2_up, base2_down, mod_depth,\
gc_rate, mismatch_density, shorest_mismatch_distance)

                #sequencing error
                (m_base1, m_base1_depth, m_base1_depth_forward,\
                        m_base1_depth_reverse, m_base2, m_base2_depth,\
                        m_base2_depth_forward, m_base2_depth_reverse) =\
                        baseSequencingAnalysis(m_base1, m_base1_depth,\
                        m_base1_depth_forward, m_base1_depth_reverse,\
                        m_base1_baq_avg, m_base1_maq_avg,base1_up, base1_down,\
                        m_base2, m_base2_depth, m_base2_depth_forward,\
                        m_base2_depth_reverse, m_base2_baq_avg,\
                        m_base2_maq_avg, base2_up, base2_down, error_snp_info)

                if (m_base1 != 'N' and m_base1 != m_reference) or\
                        (m_base2 != 'N' and m_base2 != m_reference):
                    (m_base1, m_base1_depth, m_base1_depth_forward,\
                            m_base1_depth_reverse, m_base2, m_base2_depth,\
                            m_base2_depth_forward, m_base2_depth_reverse) =\
                            baseChange(m_base1, m_base1_depth,\
                            m_base1_depth_forward, m_base1_depth_reverse,\
                            m_base2, m_base2_depth, m_base2_depth_forward,\
                            m_base2_depth_reverse, m_reference)

                    detect_depth = m_base1_depth + m_base2_depth

                    snp_rate = float(m_base2_depth) / float(detect_depth)

                    if snp_rate >= min_rate and m_base2_depth >= 2 and\
                            detect_depth <= max_depth:
                        if m_coordinate in local_indel_dic:
                            print >>sys.stderr, 'LOCALINDEL', error_snp_info
                        else:
                            if snp_rate < 0.4 and (m_base2_depth_forward < 2\
                                    or m_base2_depth_reverse < 2):
                                print >>sys.stderr, 'SNPRATE', error_snp_info
                            else:
                                all_base_info = [chrom, m_coordinate,\
                                        m_reference, m_base1, m_base1_depth,\
                                        m_base1_depth_forward,\
                                        m_base1_depth_reverse, m_base2,\
                                        m_base2_depth, m_base2_depth_forward,\
                                        m_base2_depth_reverse, snp_rate,\
                                        mismatch_density]
                                snp_for_z.append(all_base_info)

                    else:
                        print >>sys.stderr, 'DEPTH', error_snp_info
                else:
                    print >>sys.stderr, 'NOTSNP', error_snp_info

            #remove left most base info
            matchs = matchs[1:]
            coordinates = coordinates[1:]
            sequences = sequences[1:]
            mismatch_distribution = mismatch_distribution[1:]
            depths = depths[1:]

    return snp_for_z


def getWGD(wgd):
    wgd_dic = {}
    for line in open(wgd):
        line = line.rstrip()
        (chrom, coordinate, wgdup) = line.split()
        coordinate = int(coordinate)
        wgdup = int(wgdup)
        wgd_dic[coordinate] = wgdup

    return wgd_dic


def detectSNPDistance(snps):
    snp_distance_dic = {}
    for i in range(1, len(snps)):
        coordinate1 = snps[i - 1][1]
        coordinate2 = snps[i][1]
        distance = abs(coordinate2 - coordinate1) - 1
        if distance < 0:
            print >>sys.stderr, "ERROR: SNP DISTANCE < 0",\
                    coordinate1, coordinate2
            sys.exit(1)
        if coordinate1 in snp_distance_dic:
            if snp_distance_dic[coordinate1] > distance:
                snp_distance_dic[coordinate1] = distance
        else:
            snp_distance_dic[coordinate1] = distance

        if coordinate2 in snp_distance_dic:
            if snp_distance_dic[coordinate2] > distance:
                snp_distance_dic[coordinate2] = distance
        else:
            snp_distance_dic[coordinate2] = distance
    return snp_distance_dic


def joinSNPInfo(snp):
    (chrom, coordinate, reference, base1, base1_depth, base1_depth_forward,\
            base1_depth_reverse, base2, base2_depth, base2_depth_forward,\
            base2_depth_reverse, snp_rate, mismatch_density) = snp

    snp_dic = {"AT" : "W", "TA" : "W",
            "CG" : "S", "GC" : "S",
            "TG" : "K", "GT" : "K",
            "AC" : "M", "CA" : "M",
            "CT" : "Y", "TC" : "Y",
            "AG" : "R", "GA" : "R",
            'A' : 'A', 'T' : 'T', 'G' : 'G', 'C' : 'C'}

    depth = base1_depth + base2_depth
    if snp_rate > 0.9 and base1_depth <= 2:
        snp_rate = 1

    snp_base = ''
    if snp_rate == 1:
        snp_base = base2
    else:
        tem_base = "".join([base1, base2])
        if tem_base in snp_dic:
            snp_base = snp_dic[tem_base]
        else:
            print >>sys.stderr, 'JOINERROR', chrom, coordinate, reference,\
                    base1, base1_depth, base1_depth_forward,\
                    base1_depth_reverse, base2, base2_depth,\
                    base2_depth_forward, base2_depth_reverse, snp_rate,\
                    mismatch_density
            sys.exit(1)

    base1_depth_info = ":".join([str(base1_depth),\
            ",".join([str(base1_depth_forward), str(base1_depth_reverse)])])
    base2_depth_info = ":".join([str(base2_depth),\
            ','.join([str(base2_depth_forward), str(base2_depth_reverse)])])

    return "%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%.2f\t%d" % \
            (chrom, coordinate, reference, snp_base, depth, base1,\
            base1_depth_info, base2, base2_depth_info, snp_rate,\
            mismatch_density)


def snpAnalysisZ(pileup, wgd, min_baq, min_maq, min_rate, read_len, max_depth, output):
    read_len = int(read_len)
    max_depth = int(max_depth)
    min_baq = int(min_baq)
    min_maq = int(min_maq)
    min_rate = float(min_rate)

    snps = baseAnalysisY(pileup, min_baq, min_maq, min_rate, read_len, max_depth)
    snp_distance_dic = detectSNPDistance(snps)
    wgd_dic = {}
    if wgd != "-":
        wgd_dic = getWGD(wgd)

    OUTPUT = open(output, 'w')
    OUTPUT.write("#chrom\tcoordinate\treference\tSNP\tdepth\tbase1\t\
base1D:base1FD,base1RD\tbase2\tbase2D:base2FD,base2RD\trate\tmisDensity\t\
snpDistance\tWGD\n")
    for snp in snps:
        snp_info = joinSNPInfo(snp)
        coordinate = snp[1]
        coordinate_wgd = 1
        if wgd != "-":
            coordinate_wgd = 0
        elif coordinate in wgd_dic:
            coordinate_wgd = wgd_dic[coordinate]
        snp_distance = snp_distance_dic[coordinate]
        OUTPUT.write("%s\t%d\t%d\n" % (snp_info, snp_distance, coordinate_wgd))
    OUTPUT.close()


if len(sys.argv) < 2:
    print "snp.py <cns> <wgd> <min_baq> <min_maq> <min_rate> <read_len>\
 <max_depth> <snp_output>"
    sys.exit(0)
else:
    snpAnalysisZ(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],\
            sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
