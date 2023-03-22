#!/usr/bin/env python

import time, re
import sys, os, os.path
from fnmatch import fnmatch
from collections import defaultdict
from pick_up_remain import *
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

__author__ = "Xiongbin Kang"


usage = """%prog python HERO.py -r <high quality reads > -lc <long reads/contigs>  -o <output file>  -t <threads (int)> -s <split number (int)>

 This program aim to utilize high quality short or long reads to correct sequencing error for long reads or contigs. The input format of long read or contigs suggest fasta.

"""

def main():
    parser = ArgumentParser(prog='python HERO.py', description=usage, epilog='Built by: {}'.
                            format(__author__), formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', dest = 'short_reads', required=True, type=str, help='The input file of high quality short reads sequence. you also can use high quality long reads as input, but you need to activate the parameter -hlong_reads.')
    parser.add_argument('-lc', dest = 'long_reads', required=True, type=str, help='The input file of long reads/contigs (fasta) which will be polished')
    parser.add_argument('-o', dest = 'out_con', required=False, type=str, default = "polished.fa", help='The output file of corrected long reads/contigs')
    parser.add_argument('-i', dest = 'iteration', required=False, type=int, default = 3, help='The number for iterating builds overlap graph to correct target reads/contigs')
    parser.add_argument('-pl', dest = 'platform', required=False, type=str, default = "pb", help='Sequencing platform for long reads: pb/ont')
    parser.add_argument('-t', dest = 'threads', required=False, type=int, default = 30, help='the number of threads')
    parser.add_argument('-s', dest = 'split_nu', required=False, type=int, default = 10, help='The time for splitting long reads/contigs file to accelerate speed. When the high quality reads file more than 3G, we suggest split long read or contigs excess 30 times.')
    parser.add_argument('-m', dest = 'mc', required=False, type=int, default = 5, help='The mini coverage to define a SNP')
    parser.add_argument('-hlong_reads', dest = 'hlong_reads', action='store_true', help='When you utilize high quality long reads to correct long reads / contigs, please activate this parameter.')

    args = parser.parse_args()

    if not (args.short_reads and args.long_reads):
        print("Please enter long reads or contigs which will be corrected and the high quality reads that will be used to correct them.")
        parser.print_help()
        sys.exit()

    bin = os.path.split(os.path.realpath(__file__))[0]
    con = args.long_reads
    reads = args.short_reads
    out_con = args.out_con
    threads = args.threads
    mc = args.mc
    hlong_reads = args.hlong_reads
    split_nu = args.split_nu
    platform = "pb" if args.platform == "pb" else "ont"
    folder = os.getcwd()
    
    k = 0
    while(k < args.iteration):
        k += 1
        out_con2 = folder + "/" + str(k) + "_" + out_con
        
        correct(con, reads, out_con2, folder, mc, split_nu, threads, hlong_reads, bin, platform)
        con = out_con2

    out_con2 = folder + "/" + "final_" + out_con
    
    if hlong_reads:
        
        execute("mv %s %s;" %(con, out_con2))
        
    else:

        hlong_reads = True
        reads = os.path.basename(con)
        correct(con, reads, out_con2, folder, mc, split_nu, threads, hlong_reads, bin, platform)


def correct(con, reads, out_con2, folder, mc, split_nu, threads, hlong_reads, bin, platform):

        sub = "sub"+str(time.time())[-3:]
        folder2 = folder + "/" + sub + "_intermedia"
        execute("mkdir -p %s" %folder2)

        fx = fq_or_fa(con)
        con = fq2fa(con, fx, folder2)
        fx = fq_or_fa(con)
        th2 = 3 + int(threads/split_nu)

        nu = int(os.popen("wc -l %s" %con ).readline().split()[0])
        nu_sub = int(nu/(8*split_nu)+1)*8

        split_line = "cd %s; split %s -l %s -d -a 3 %s;" %(folder2, con, nu_sub, sub)
        execute(split_line)

        # run polisd in small group
        sub_overlap = [name for name in os.listdir(folder2) if fnmatch(name, sub+"*")]
        cmd_rename = ["mv %s %s.%s" %(i, i, fx) for i in sub_overlap]
        cmd_rename2 = ";".join(cmd_rename)
        execute2(("cd %s;" + cmd_rename2)%folder2)
    
        sub_overlap2 = [name for name in os.listdir(folder2) if fnmatch(name, sub+"*")]
    
        if hlong_reads:

            cmd = ["python %s/pre_polish.py -pl %s -r %s/%s -c %s/%s -mc %s -t %s -long_reads " %(bin, platform, folder, reads, folder2, i, mc, th2) for i in sub_overlap2]
        
        else:
            cmd = ["python %s/pre_polish.py -r %s/%s -c %s/%s -mc %s -t %s " %(bin, folder, reads, folder2, i, mc, th2) for i in sub_overlap2]

        cmd2 = "\n".join(cmd)
        
        with open('%s/%s_overlap.sh' %(folder2, sub), 'w') as sh:
            sh.write(cmd2)
            sh.write("\n")
        
        threads2 = 10 if threads > 10 else threads
        run_cmd = "cd %s; cat %s_overlap.sh | xargs -i -P %s bash -c \"{}\";" %(folder2, sub, threads2)
        execute(run_cmd)
        execute("cd %s; cat %s*_polished.fa > %s;" %(folder2, sub, out_con2))
    
        # pick up the reads that weren't corrected.
        remain_fa = folder2 + "/%s_all_remain.fa" %(sub)
        file_model = fq_or_fa2(con)
        pick_remain(con, out_con2, remain_fa, file_model)
        execute("cd %s; cat %s >> %s; cd %s; rm -r %s*;" %(folder2, remain_fa, out_con2, folder, sub))
  

def fq2fa(fq, model, folder2):

    """Filter out non actg base"""
    i = 0
    new_out_fa = folder2+"/pre.fa"

    with open(new_out_fa, 'w') as out1:
    
        if model == "fastq":
            with open(fq,'r') as file:
                for num, line in enumerate(file):
                    if num % 4 == 1:
                        newline = re.sub(r'[^ATGCN\n]', "N", line.upper())
                        out1.write(newline)
                    elif num % 4 == 0 and line.startswith('@'):
                        
                        new_id1 = line.strip().split(" ")[0]
                        new_id = ">" + new_id1[1:] + "\n"
                        out1.write(new_id)
        else:
            with open(fq,'r') as file:
                for num, line in enumerate(file):
                    if num % 2 == 1:
                        newline = re.sub(r'[^ATGCN\n]', "N", line.upper())
                        out1.write(newline)
                    else:
                        new_id1 = line.strip().split(" ")[0]
                        new_id = new_id1 + "\n"
                        out1.write(new_id)           
        
    return(new_out_fa)

def execute(cmd):
    te = os.system(cmd + " 2> log_error.txt")
    if te:
        with open("log_error.txt","r") as file:
            print("Don't execute the command: %s " %cmd)
            print(file.read())
    else:
        print("successfully execute: %s" %cmd)

def execute2(cmd):
    te = os.system(cmd + " 2>log_error.txt")
    if te:
        with open("log_error.txt","r") as file:
            print("Don't execute the command: %s " %cmd)
            print(file.read())
        
def fq_or_fa(file):
    # with open(file) as fr:
        # s = fr.readline()[0]
    s = os.popen("less {}|head -1".format(file)).read()[0]
    mode = ''
    if s == '>':
        mode = 'fasta'
    elif s == '@':
        mode = 'fastq'
    else:
        raise Exception("invalid input file, must be FASTA/FASTQ format.", file)
    return mode



if __name__ == '__main__':
        sys.exit(main())
