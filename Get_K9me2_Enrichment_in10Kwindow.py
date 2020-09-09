#!/usr/bin/python
# -*-coding:utf-8-*-
# @author xialinjian

import os
from multiprocessing import Pool
from scipy import stats
import pandas as pd
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

Bam_dir ='' 
chrom_size=''
Result_dir=''
Input_reads_limit = 100
window_size=10000

def get_prefix():
    line=[i.split('.')[0] for i in os.listdir("%s/IP"%(Bam_dir)) if i.endswith('.bam')]
    return line

def get_bam(prefix):
    IP,Input = "%s/IP/%s.bam"%(Bam_dir,prefix),"%s/Input/%s.bam"%(Bam_dir,prefix)
    return IP,Input

def make_window():
    window="%s/window_%d.bed"%(Result_dir,window_size)
    os.system("bedtools makewindows -g %s -w %d|grep -v '_' > %s"%(chrom_size,window_size,window))
    return window

def get_split(window):
    tmp = "%s/chr_split" %(Result_dir)
    if not os.path.exists(tmp):
        os.mkdir(tmp)
    with open("%s"%(window),'r') as f:
        for eachline in f:
            if eachline:
                chrom = eachline.strip("\n").split('\t',1)[0]
                if chrom.startswith('chr'):
                    with open("%s/%s.txt" %(tmp,chrom),'a+') as out:
                        out.write(eachline)
                else:
                    continue
            else:
                break

def get_line():
    line=[]
    for i in os.listdir("%s/chr_split"%(Result_dir)):
        if i.endswith('.txt'):
            line.append(i)
    return line

def get_padjust(pvalue_list):
    stats = importr('stats')
    p_adjust = stats.p_adjust(FloatVector(pvalue_list), method = 'BH')
    return p_adjust

def get_coverage(chr,prefix):
    IP,Input=get_bam(prefix)
    bed = "%s/chr_split/%s" % (Result_dir, chr)
    result_dir = "%s/%s" % (Result_dir, prefix)
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
    os.system("bedtools multicov -bams %s %s -bed %s > %s/%s.coverage" %(IP,Input,bed,result_dir,chr))

def get_total_reads(file):
    with open("%s"%(file),'r') as f:
        line=f.readlines()[0]
        number = int(line.rstrip().split()[0])
    return number

def get_IP_Input_sum(prefix):
    IP, Input = get_bam(prefix)
    IP_log,Input_log = "%s/IP/%s.log"%(Bam_dir,prefix),"%s/Input/%s.log"%(Bam_dir,prefix)
    os.system("samtools flagstat %s > %s"%(IP,IP_log))
    os.system("samtools flagstat %s > %s"%(Input,Input_log))
    IP_sum,Input_sum = get_total_reads(IP_log),get_total_reads(Input_log)
    return IP_sum,Input_sum


def get_fisher(chr,prefix):
    IP_sum,Input_sum = get_IP_Input_sum(prefix)
    coverage_file = "%s/%s/%s.coverage" % (Result_dir, prefix, chr)
    get_coverage(chr,prefix)
    fisher_txt = "%s/%s/%s.result"%(Result_dir,prefix,chr)
    with open("%s" %(fisher_txt),'w') as fw:
        with open("%s" %(coverage_file),'r') as f:
            a =f.readlines()
        for i in a:
            j = i.rstrip().split()
            stat,p = stats.fisher_exact([[int(j[3]),int(j[4])],[IP_sum,Input_sum]])
            if int(j[3]) != 0 and IP_sum != 0 and Input_sum != 0:
                FC=str((float(j[3])/IP_sum)/(float(j[4])/Input_sum))
            else:
                FC="NA"
            fw.write("%s\t%d\t%d\t%s\t%f\n"%(i.rstrip(),IP_sum,Input_sum,FC,p))


def main():
    window=make_window()
    get_split(window)
    sample_line = get_prefix()
    for sample in sample_line:
        IP_sum, Input_sum = get_IP_Input_sum(sample)
        chr_line=get_line()
        pool=Pool()
        for chr in chr_line:
            pool.apply_async(get_fisher, (chr, sample))
        pool.close()
        pool.join()

if __name__=='__main__':
    main()



