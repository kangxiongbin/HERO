#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
from collections import defaultdict
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

__author__ = "Xiongbin Kang"


usage = """%prog python pre_polish.py -r <reads >  -c <contigs> -mc <int/mini coverage> -t <int/threads>

 This program is used reads to polish contigs (refine racon)

"""

def main():
    parser = ArgumentParser(description = usage)
    parser.add_argument('-r', dest = 'reads', required=True, type=str, help='sequence')
    parser.add_argument('-c', dest = 'con', required=True, type=str, help='contigs')
    parser.add_argument('-mc', dest = 'mc', required=False, type=int, default = 5, help='the mini coverage of one mutation (snp)')
    parser.add_argument('-t', dest = 'threads', required=False, type=int, default = 10, help='the number of threads')
    parser.add_argument('-pl', dest = 'platform', required=False, type=str, default = "pb", help='sequencing platform: pb/ont')
    parser.add_argument('-min_len', dest = 'min_len', required=False, type=int, default = 3000, help='the min length of overlap region')
    parser.add_argument('-min_mis', dest = 'min_mis', required=False, type=float, default = 0.02, help='the min mismatch in the overlap region')
    parser.add_argument('-long_reads', dest = 'long_reads', action='store_true', help='Long reads activte the parameter')
    parser.add_argument('-thre', dest = 'threshold', required=False, default = 0.0005, type = int, help='allow the number of mutations in the overlap region between two reads. it just work for long reads')
    
    args = parser.parse_args()

    if not (args.reads and args.con):
        print("Please input reads and contigs")
        parser.print_help()
        sys.exit()

    bin = os.path.split(os.path.realpath(__file__))[0]
    con = args.con
    reads = args.reads
    mc = args.mc
    threads = args.threads
    min_len = args.min_len
    min_mis = args.min_mis
    threshold = args.threshold
    platform = "pb" if args.platform == "pb" else "ont"

    
#    folder = os.getcwd()
    overlap_long = "%s_tmp_long_overlap.paf" %(con)
    overlap = "%s_tmp_overlap.paf" %(con)
    overlap_filter = "%s_tmp_overlap2.paf" %(con)
    
    overlap_sort = "%s_tmp_sorted_overlap.paf" %(con)
	
    overlap2 = "%s_overlap.paf" %(con)

    if args.long_reads:
        execute("minimap2 -t %s -L --eqx -cx ava-%s -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25 %s %s 2> /dev/null 1> %s" %(threads, platform, con, reads, overlap_long))
        execute("minimap2 -t %s -L --eqx -cx ava-%s -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25 %s %s 2> /dev/null 1>> %s" %(threads, platform, reads, con, overlap_long))
        
        paf_ovlp1 = os.popen("cat %s" %overlap_long)
        execute("rm %s"%overlap_long)
        with open (overlap, "w") as olvp_out:
            dedup = {}
            for line in paf_ovlp1:
                len_sp = line.split('\t')
                
                if int(len_sp[9]) < min_len:
                    continue
                
                mis = int(len_sp[12].split(':')[2]) / int(len_sp[9])
                if mis > min_mis:
                    continue
                
                k1 = ':'.join(sorted([len_sp[0], len_sp[5]]))
                if k1 in dedup:
                    continue 
                else:
                    dedup[k1] = 1

                olvp_out.write(line)

    else:
    
        execute("minimap2 -t %s  -L --eqx -c --sr -D --no-long-join -k 21 -w 11 -s 60 -m 30 -n 2 -A 4 -B 2 -N 2 --end-bonus=100 %s %s 2> /dev/null 1> %s" %(threads, con, reads, overlap))

    count = 0
    
    size = os.path.getsize(overlap)

#    if size > 3*1024*1024*1024:
    
#        overlap = filter_file(overlap, overlap_filter)

    while True:
        execute("sort -nk7 -k8 -k9 -k5 %s > %s;" % (overlap, overlap_sort))   
        
        if os.path.exists(overlap_sort) and os.path.getsize(overlap_sort) > 0:
            execute("rm %s; " %(overlap))
            break
            
        count += 1
        if count > 5:
            break
	
    count = 0
    while True:
    
        paf_file2 = os.popen("cat %s" %overlap_sort)

        if args.long_reads:
    
            (snp, map_po, start_po) = prpare_mutation2(paf_file2)
        
        else:
    
            (snp, map_po, start_po) = prpare_mutation(paf_file2)
 
        start_po_sorted = defaultdict(list) # read start and end
        
        for k in start_po.keys():
            start_po_sorted[k] = sorted(start_po[k], key = lambda x: (x[0], x[1]))
	
        mutation = mutation_re(snp, start_po_sorted, map_po, mc = mc)
	
        mutation2 = defaultdict(list)

        paf_file2 = os.popen("cat %s" %overlap_sort)

        with open (overlap2, "w") as out:
            if args.long_reads:
                for line1 in paf_file2:
                    [qseqid, qlen, qstart, qend, qori, sseqid, slen, sstart, send, matchcount, length] = line1.split('\t')[:11]
                    fkey = ":".join(sorted([qseqid, sseqid]))
                    if fkey in mutation:
                        fre = mutation[fkey]/int(matchcount)
                        if fre > threshold:
                            continue
                        
                    out.write(line1)
        
            else:
                for line1 in paf_file2:
                    line = line1.split("\t")
                    fkey = ":".join(sorted([line[0], line[5]]))

                    if fkey in mutation:
                        continue
                    else:
                        out.write(line1)
                        
        if os.path.exists(overlap2) and os.path.getsize(overlap2) > 0:
    
            execute("rm %s;"%overlap_sort)
            break

        count += 1
        if count > 6:
            break

    polished_fa = "%s_polished.fa" %(con)
    remain_fa = "%s_remain.fa" %(con)
    
    logfile = "%s_logfile" %(con)
    threads = 5
    execute("racon --no-trimming -u -t {} {} {} {} >{} 2>{};".format(threads, reads, overlap2, con, polished_fa, logfile))
#    execute("racon -t {} {} {} {} >{} 2>{};".format(threads, reads, overlap2, con, polished_fa, logfile))

    # if the output of racon exists, remove Intermediate redundant file 
    if os.path.exists(polished_fa) and os.path.getsize(polished_fa):
    
        execute("rm %s %s"%(overlap2, logfile))
            
    if args.long_reads:
    
        pick_remain(polished_fa, con, remain_fa)
        
        execute("cat %s >> %s" %(remain_fa, polished_fa))


def filter_file(in_file, out_file):

    name_counts = {}
    with open(in_file) as f_in, open(out_file, 'w') as f_out:
        for line in f_in:
            name = line.split()[0]

            if name not in name_counts:
                name_counts[name] = 0

            name_counts[name] += 1

            if name_counts[name] <= 5:
                f_out.write(line)

        execute("rm %s; " %(in_file))
        return out_file
 
def pick_remain(polished_fa, con, remain_fa):

    d = defaultdict(list)
    
    with open(polished_fa) as f:
        i = 0
        for line in f:
            if i % 2 == 0 :
                line.replace('\n','')
                key = line.split( )[0]
                d[key] = 1

            i = i + 1
        
    with open(remain_fa,'a') as fw:
        with open(con) as q:
    
            i = 0
            for line2 in q:
                if i % 2 == 0 :
                    line2.replace('\n','')
                    key2 = line2.split( )[0]
                    if key2 in d:
                        p1 = 0
                    else:
                        p1 = 1

                if p1 > 0:
                    fw.write(line2)
                i = i + 1

def plus_dict(key, dict):
    
    dict[key] = dict[key] + 1 if key in dict else 1
    
    return(dict)

def plus_value(key, value, dict):
    
    if key in dict:
        dict[key].append(value)
                    
    else:
        dict[key] = [value]
    
    return(dict)

def plus_value2(key, value1, value2, dict):
    
    if key in dict:
        dict[key].append([value1, value2])
                    
    else:
        dict[key] = [[value1, value2]]
    
    return(dict)
    
def execute(cmd):
    te = os.system(cmd + " 2>out.txt")
    if te:
        with open("out.txt","r") as file:
            print("Don't execute the command: %s " %cmd)
            print(file.read())
        
def prpare_mutation(paf_file):

        snp = defaultdict(list) # snp position and times
        map_po = defaultdict(list) # snp position and the read that support the snp
        start_po = defaultdict(list) # the start and end position of match reads
        deduplication = defaultdict(list) # delete the duplication of read combination

        ref_len = {}

        for line in paf_file:
            len_sp = line.split('\t')
            
            if len(len_sp) < 6:
                continue

            qseqid = len_sp[0]
            flag = len_sp[4]
            refid = len_sp[5]
            refposi = len_sp[7]
            refposi2 = len_sp[8]
            cigar1 = len_sp[-1]

            if qseqid == refid:
                continue
                
            if cigar1 == "*":
                continue
                
            cigar = cigar1.strip()
            ty =  re.split("\d+", cigar) # the array contains all type of mutation
            po1 = re.split("\D", cigar) # # the array contains the number of mutation (how many based in a mutation)

            ty.pop(0)
            po = po1[5:-1]

            pos2 = int(refposi)

            start_po = plus_value2(refid, int(refposi), int(refposi2), start_po)
                
            for t in range(len(ty)):
                
                if ty[t] == "=":
                    pos2 += int(po[t])
 
                elif ty[t] == "D":
                    pos2 += int(po[t])
                    
                elif ty[t] == "X":

                    pos2 += int(po[t])

                    key2 = ":".join([refid, str(pos2)])

                    snp = plus_dict(key2, snp)

                    map_po = plus_value(key2, qseqid, map_po)
                    
                    
        return(snp, map_po, start_po)

def prpare_mutation2(paf_file):

        snp = defaultdict(list) # snp position and times
        map_po = defaultdict(list) # snp position and the read that support the snp
        start_po = defaultdict(list) # the start and end position of match reads
        deduplication = defaultdict(list) # delete the duplication of read combination

        ref_len = {}

        for line in paf_file:
            len_sp = line.split('\t')
            
            if len(len_sp) < 6:
                continue
                
            qseqid = len_sp[0]
            qseqlen = len_sp[1]
            qseqposi = len_sp[2]
            qseqposi2 = len_sp[3]
            flag = len_sp[4]
            refid = len_sp[5]
            refposi = len_sp[7]
            refposi2 = len_sp[8]
            cigar1 = len_sp[-1]
            cigar = cigar1.strip()
            
            if qseqid == refid:
                continue
                
            if cigar == "*":
                continue
                
            k1 = ':'.join(sorted([qseqid, refid]))
            
            if k1 in deduplication:
                continue        
            else:
                deduplication[k1] = 1
            
            ty =  re.split("\d+", cigar) # the array contains all type of mutation
            po1 = re.split("\D", cigar) # # the array contains the number of mutation (how many based in a mutation)

            ty.pop(0)
            po = po1[5:-1]
            
            pos1 = int(qseqposi) if flag == "+" else int(qseqlen) - int(qseqposi2)
            
            pos2 = int(refposi)

            start_po = plus_value2(refid, int(refposi), int(refposi2), start_po)
            start_po = plus_value2(qseqid, int(qseqposi), int(qseqposi2), start_po)

            for t in range(len(ty)):
                
                if ty[t] == "=":
                    pos1 += int(po[t])
                    pos2 += int(po[t])

                elif ty[t] == "I":
                    pos1 += int(po[t])
                    
                elif ty[t] == "D":
                    pos2 += int(po[t])
                    
                elif ty[t] == "X":

                    pos1 += int(po[t])
                    pos2 += int(po[t])

                    key1 = ":".join([qseqid, str(pos1)]) if flag == "+" else ":".join([qseqid, str(int(qseqlen) - pos1 + 1)])
                    key2 = ":".join([refid, str(pos2)])

                    snp = plus_dict(key1, snp)
                    snp = plus_dict(key2, snp)
                    
                    map_po = plus_value(key1, refid, map_po)
                    map_po = plus_value(key2, qseqid, map_po)                    

        return(snp, map_po, start_po)

def mutation_re(snp, start_po, map_po, mc = 5):

    mutation = {}

    for k,v in snp.items():
        
        if v < mc: # at least three reads support SNP
                continue

        (read, position) = k.split(":")
        position = int(position)
        
        con = 0
        for p in start_po[read]:
            
            if position > p[1]:
                continue

            elif position > p[0] and position < p[1] :
                con += 1

            elif position < p[0]:
                break
				
        c = con - v 

        if c < mc: # how many reads support the query reads base type
             continue
            
        for ke in map_po[k]:
                
              fkey = ':'.join(sorted([ke, read]))
                
              mutation = plus_dict(fkey, mutation)
 
    return(mutation)

if __name__ == '__main__':
        sys.exit(main())
