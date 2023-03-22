#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
from collections import defaultdict
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

__author__ = "Xiongbin Kang"


usage = """%prog python pick_up_remain.py -p <previou file >  -c <current file> -r <remain reads> -t <int/threads>

 This program is used to pick up remianing reads

"""

def main():
    parser = ArgumentParser(description = usage)
    parser.add_argument('-p', dest = 'pre_fa', required=True, type=str, help='previou sequence file')
    parser.add_argument('-c', dest = 'cur_fa', required=True, type=str, help='current sequence file')
    parser.add_argument('-r', dest = 'remain_fa', required=False, type=str, default = "remain.fa", help='remaining sequence in the file as output')
    args = parser.parse_args()

    if not (args.pre_fa and args.cur_fa):
        print("Please input reads file")
        parser.print_help()
        sys.exit()
        
    pre_fa = args.pre_fa
    cur_fa = args.cur_fa
    remain_fa = args.remain_fa
    file_model = fq_or_fa2(cur_fa)
    pick_remain(pre_fa, cur_fa, remain_fa, file_model)
    
def pick_remain(pre_fa, cur_fa, remain_fa, file_model):

    d = defaultdict(list)
    t = int(file_model)
    
    with open(cur_fa) as f:
        i = 0
        for line in f:
            if i % t == 0 :
                line.replace('\n','')
                key = line.split( )[0]
                d[key] = 1

            i = i + 1
        
    with open(remain_fa,'w') as fw:
        with open(pre_fa) as q:
    
            j = 0
            for line2 in q:
                if j % t == 0 :
                    line2.replace('\n','')
                    key2 = line2.split( )[0]
                    
                    if key2 in d:
                        p1 = 0
                    else:
                        p1 = 1

                if p1 > 0:
                    
                    if(j % t == 1):
                        fw.write(line2.upper())
                    else:
                        fw.write(line2)
                    
                j = j + 1
                
    p = round(i/j, 2)
     
    if(p < 0.85):
       print("Only %s reads were corrected. If the coverage more than 5, it means that some sub-threads did not complete. We recommend increase the number of split and rerun this program. If you are correcting contigs, you can ignore this warning." %(p))
        
    else:
       print("%s reads were corrected in this cycle." %(p))


def fq_or_fa2(file):
    # with open(file) as fr:
    # s = fr.readline()[0]
    s = os.popen("less {}|head -1".format(file)).read()[0]
    mode = ''
    if s == '>':
        mode = 2
    elif s == '@':
        mode = 4
    else:
        raise Exception("invalid input file, must be FASTA/FASTQ format.", file)
    return mode


if __name__ == '__main__':
        sys.exit(main())
        