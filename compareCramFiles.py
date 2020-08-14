###################################
##  compare crams parallel
##  read by read
##  write reads aligend differently
##  to file
###################################

import pysam
import hashlib
import time
from multiprocessing import Pool
import sys
from itertools import zip_longest
import os
import argparse



class CustomFormatter(argparse.RawDescriptionHelpFormatter):
    pass

epilog="""
Program that compares two cram files and checks whether or not their content is the same
the script accepts 6 command line arguments:
-i1, -i2, the two input cram files respectively for comparison
-t file type for the input files <s, b, c> for sam, bam and cram respectively.
-r reference file directory used for the read maps
-n number of threads used for multiprocessing
-o output directory to store the discrepencies between the two files in comparison if any
"""

parser = argparse.ArgumentParser(description="CRAM files comparison tool",epilog=epilog, formatter_class=CustomFormatter)
parser.add_argument("-i1", "--cram_f1",
                    action="store",
                    dest="cram_f1",
                    required=True,
                    help=("directory for the first cram file"))
parser.add_argument("-i2", "--cram_f2",
                    action="store",
                    dest="cram_f2",
                    required=True,
                    help=("directory for the second cram file"))
parser.add_argument("-t", "--file_type",
                    action="store",
                    dest="file_type",
                    required=True,
                    help=("type of the files in comparison"))
parser.add_argument("-r", "--reference_fname",
                    action="store",
                    dest="reference_fname",
                    required=True,
                    help=("directory for the reference genome used to map the reads"))
parser.add_argument("-n", "--n_processes",
                    action="store",
                    dest="n_processes",
                    required=True,
                    help=("number of threads used to run the script"))
parser.add_argument("-o", "--out_dir",
                    action="store",
                    dest="out_dir",
                    required=False,
                    help=("output directory to store the files"))

args = parser.parse_args()

if not os.path.isfile(args.cram_f1) or not os.path.isfile(args.cram_f2):
    raise RuntimeError("Nonexistant or invalid input directory '{}'".format(args.cram_f1))

cram_f1 = args.cram_f1
cram_f2 = args.cram_f2
file_type = args.file_type
reference_fname = args.reference_fname
n_processes = int(args.n_processes)
out_dir = args.out_dir

ftype_dic = {'-s':'r', '-b':'rb', '-c':'rc', 'c':'rc', 's':'r', 'b':'rb'}

file_type = ftype_dic[file_type]

def getMD5digest(cram_f, byte_size = 8191):
    with open(cram_f, "rb") as f:
        file_hash = hashlib.md5()
        while True:
            chunk = f.read(byte_size)
            if not chunk:
                break
            file_hash.update(chunk)
    return file_hash.hexdigest()

def cramSectionReadByRead(cram_section1, cram_section2):
    seq_hash_set1 = set()
    seq_hash_set2 = set()
    for read1, read2 in zip_longest(cram_section1, cram_section2):
        if (read1 != None): #EOF
            qname1, flag1, ref_start1 = read1.query_name, str(read1.flag), str(read1.reference_start)
            if read1.cigarstring == None:
                cigar1 = '-1'
            else:
                cigar1 = read1.cigarstring
            read_lst1 = [qname1, flag1, ref_start1, cigar1]
            seq_str1 = ''.join(read_lst1)
            seq_hash1 = hashlib.sha1(str.encode(seq_str1)).hexdigest()

        if(read2 != None):
            qname2, flag2, ref_start2 = read2.query_name, str(read2.flag), str(read2.reference_start)
            if read2.cigarstring == None:
                cigar2 = '-1'
            else:
                cigar2 = read2.cigarstring
            read_lst2 = [qname2, flag2, ref_start2, cigar2]
            seq_str2 = ''.join(read_lst2)
            seq_hash2 = hashlib.sha1(str.encode(seq_str2)).hexdigest()

        if seq_hash1 == seq_hash2:
            continue
        if seq_hash1 in seq_hash_set2:
            seq_hash_set2.remove(seq_hash1)
            
        elif seq_hash1 not in seq_hash_set2:
            seq_hash_set1.add(seq_hash1)
            
            
        if seq_hash2 in seq_hash_set1:
            seq_hash_set1.remove(seq_hash2)
        
        elif seq_hash2 not in seq_hash_set1:
            seq_hash_set2.add(seq_hash2)

    return len(seq_hash_set1), len(seq_hash_set2)

def cram2hash_compareChromosomes(chromosome, cram_f1, cram_f2, file_type, reference_fname):
    """
    if you want to fetch chromosome, you can always open a new file everytime here
    but you need to take care of the unmapped reads.
    check out splitting the cram file into chunks using samtools
    """
    cramfile1 = pysam.AlignmentFile(cram_f1, file_type, reference_filename = reference_fname)
    cram_section1 = cramfile1.fetch(chromosome)

    cramfile2 = pysam.AlignmentFile(cram_f2, file_type, reference_filename = reference_fname)
    cram_section2 = cramfile2.fetch(chromosome)

    nreads1, nreads2 = cramSectionReadByRead(cram_section1, cram_section2)    
    
    return (nreads1, nreads2)

def cramSectionReadByRead_output(cram_section1, cram_section2):
    seq_hash_table1 = dict()
    seq_hash_table2 = dict()
    seqID2hash_table1 = dict()
    seqID2hash_table2 = dict()
    for read1, read2 in zip_longest(cram_section1, cram_section2):
        if (read1 != None): #EOF
            qname1, flag1, ref_start1 = read1.query_name, str(read1.flag), str(read1.reference_start)        
            if read1.cigarstring == None:
                cigar1 = '-1'
            else:
                cigar1 = read1.cigarstring
            read_lst1 = [flag1, ref_start1, cigar1]
            hash_lst1 = [qname1, flag1, ref_start1, cigar1]
            seq_str1 = ''.join(hash_lst1)
            seq_hash1 = hashlib.sha1(str.encode(seq_str1)).hexdigest()
        
        if(read2 != None):
            qname2, flag2, ref_start2 = read2.query_name, str(read2.flag), str(read2.reference_start)
            if read2.cigarstring == None:
                cigar2 = '-1'
            else:
                cigar2 = read2.cigarstring
            read_lst2 = [flag2, ref_start2, cigar2]
            hash_lst2 = [qname2, flag2, ref_start2, cigar2]
            seq_str2 = ''.join(hash_lst2)
            seq_hash2 = hashlib.sha1(str.encode(seq_str2)).hexdigest()
        
        if seq_hash1 == seq_hash2:
            continue
        else:            
            if seq_hash1 in seq_hash_table2:
                del seq_hash_table2[seq_hash1]
                if qname1 in seqID2hash_table2:
                    seqID2hash_table2[qname1].remove(seq_hash1)
                    if not seqID2hash_table2[qname1]:
                        del seqID2hash_table2[qname1]
                        
            elif seq_hash1 not in seq_hash_table2:
                seq_hash_table1[seq_hash1] = '|'.join(read_lst1)
                if (qname1 not in seqID2hash_table2) or (qname1 in seqID2hash_table2 and seq_hash1 not in seqID2hash_table2[qname1]):
                    if qname1 not in seqID2hash_table1:
                        seqID2hash_table1[qname1] = [seq_hash1]                        
                    else:
                        seqID2hash_table1[qname1].append(seq_hash1)
                        
                        
            if seq_hash2 in seq_hash_table1:
                del seq_hash_table1[seq_hash2]
                if qname2 in seqID2hash_table1:
                    seqID2hash_table1[qname2].remove(seq_hash2)
                    if not seqID2hash_table1[qname2]:
                        del seqID2hash_table1[qname2]
                 
            elif seq_hash2 not in seq_hash_table1:
                seq_hash_table2[seq_hash2] = '|'.join(read_lst2)
                if (qname2 not in seqID2hash_table1) or (qname2 in seqID2hash_table1 and seq_hash2 not in seqID2hash_table1[qname2]):
                    if qname2 not in seqID2hash_table2:
                        seqID2hash_table2[qname2] = [seq_hash2]                        
                    else:
                        seqID2hash_table2[qname2].append(seq_hash2)            

    return seq_hash_table1, seqID2hash_table1, seq_hash_table2, seqID2hash_table2
        
def cram2hash_compareChromosomes_output(chromosome, cram_f1, cram_f2, file_type, reference_fname, out_dir):
    """
    if you want to fetch chromosome, you can always open a new file everytime here
    but you need to take care of the unmapped reads.
    check out splitting the cram file into chunks using samtools
    """    
    
    cramfile1 = pysam.AlignmentFile(cram_f1, file_type, reference_filename = reference_fname)
    cram_section1 = cramfile1.fetch(chromosome)    
    
    cramfile2 = pysam.AlignmentFile(cram_f2, file_type, reference_filename = reference_fname)
    cram_section2 = cramfile2.fetch(chromosome)
    
    seq_hash_table1, seqID_hash_table1, seq_hash_table2, seqID_hash_table2 = cramSectionReadByRead_output(cram_section1, cram_section2)
    
    
    
    readsMappedDifferently = set(list(seqID_hash_table1.keys())).intersection(set(seqID_hash_table2.keys()))
    
    n_readsMappedDifferently = len(readsMappedDifferently)
    
    cram_f1_onlyMapped_reads = set(list(seqID_hash_table1.keys())).difference(readsMappedDifferently)
    cram_f2_onlyMapped_reads = set(list(seqID_hash_table2.keys())).difference(readsMappedDifferently)
    
    all_rows = list()
    for readID in readsMappedDifferently:
        cram1_hashes = seqID_hash_table1[readID]
        cram2_hashes = seqID_hash_table2[readID]
        cram1Alignments = '#'.join([seq_hash_table1[cramHash] for cramHash in cram1_hashes if cramHash in seq_hash_table1])
        cram2Alignments = '#'.join([seq_hash_table2[cramHash] for cramHash in cram2_hashes if cramHash in seq_hash_table2])
        #if cram1Alignments != '' and cram2Alignments != '': #just output the reads that are mapped differently (delete this line if you want to output all)            
        row = [readID, cram1Alignments, cram2Alignments]
        all_rows.append('\t'.join(row))   
    if all_rows:
        with open(out_dir + cram_f1.rsplit('/', 1)[-1].rsplit('.', 1)[0] + '_' + cram_f2.rsplit('/', 1)[-1].rsplit('.', 1)[0]+'_differentMappedReads.tsv', 'a') as out_f:
            out_f.write('\n'.join(all_rows)+'\n')
    
    
    
    return (len(cram_f1_onlyMapped_reads), len(cram_f2_onlyMapped_reads), n_readsMappedDifferently)

def cram2hash_parallel(cram_f1, cram_f2, n_processes, file_type = file_type, reference_fname = reference_fname, out_dir = out_dir):
    
    cramfile = pysam.AlignmentFile(cram_f1, file_type, reference_filename = reference_fname)
    chromosomes = cramfile.references
    
    chromosomeHashComparison_set = set()
    
    if out_dir == None:
        with Pool(processes=n_processes) as pool:
            chr_hash_comparison_lst = pool.starmap(cram2hash_compareChromosomes, [(chrom, cram_f1, cram_f2, file_type, reference_fname) for chrom in chromosomes])
        total_nreads1, total_nreads2 = sum([item[0] for item in chr_hash_comparison_lst]), sum([item[1] for item in chr_hash_comparison_lst])
        if (total_nreads1 == 0 and total_nreads2 == 0):
            return (False, total_nreads1, total_nreads2)
        else:
            return (True, total_nreads1, total_nreads2)
        
    elif out_dir != None:
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        
        with open(out_dir + cram_f1.rsplit('/', 1)[-1].rsplit('.', 1)[0] + '_' + cram_f2.rsplit('/', 1)[-1].rsplit('.', 1)[0]+'_differentMappedReads.tsv', 'w') as out_f:
            out_f.write('ReadID\t'+cram_f1.rsplit('/', 1)[-1].rsplit('.', 1)[0]+'_FLAG|CHROM_START|CIGAR\t'+ cram_f2.rsplit('/', 1)[-1].rsplit('.', 1)[0]+'_FLAG|CHROM_START|CIGAR\n')
            
        with Pool(processes=n_processes) as pool:
            chr_hash_comparison_lst = pool.starmap(cram2hash_compareChromosomes_output, [(chrom, cram_f1, cram_f2, file_type, reference_fname, out_dir) for chrom in chromosomes])
            
        total_nreads1, total_nreads2, total_nreads_mappedDifferently = sum([item[0] for item in chr_hash_comparison_lst]), sum([item[1] for item in chr_hash_comparison_lst]), sum([item[2] for item in chr_hash_comparison_lst])
        if (total_nreads1 == 0 and total_nreads2 == 0 and total_nreads_mappedDifferently == 0):
            return (False, total_nreads1, total_nreads2, total_nreads_mappedDifferently)
        else:
            return (True, total_nreads1, total_nreads2, total_nreads_mappedDifferently)
        
        
print('##################################################')
print('comparing the files', cram_f1, cram_f2, 'using', n_processes, 'processes')

total_time = time.time()

time_start = time.time()
cram_f1_MD5 = getMD5digest(cram_f1)
time_end = time.time() - time_start
print('Finished processing MD5 for file',  cram_f1, 'in', time_end, 'seconds')

time_start = time.time()
cram_f2_MD5 = getMD5digest(cram_f2)
time_end = time.time() - time_start
print('Finished processing MD5 for file',  cram_f2, 'in', time_end, 'seconds')

time_start = time.time()
print('Are the MD5 digests for the 2 cram files the same?: ', cram_f1_MD5 == cram_f2_MD5, time.time() - time_start, 'seconds')

time_start = time.time()
cram_chromosomeHashTable_comparisons = cram2hash_parallel(cram_f1, cram_f2, n_processes)


time_end = time.time() - time_start
print('Finished processing the cram files',  cram_f1, cram_f2, 'in', time_end, 'seconds')

print('Are the two cram files in comparison different?', cram_chromosomeHashTable_comparisons[0])
if out_dir == None:
    if cram_chromosomeHashTable_comparisons[0] == True:
        print('there are', cram_chromosomeHashTable_comparisons[1], 'reads in', cram_f1.rsplit('/', 1)[-1], 'not found mapped in', cram_f2.rsplit('/', 1)[-1])
        print('there are', cram_chromosomeHashTable_comparisons[2], 'reads in', cram_f2.rsplit('/', 1)[-1], 'not found mapped in', cram_f1.rsplit('/', 1)[-1])
        
if out_dir != None:
    if cram_chromosomeHashTable_comparisons[0] == True:
        print('there are', cram_chromosomeHashTable_comparisons[1], 'reads in', cram_f1.rsplit('/', 1)[-1], 'not found mapped in', cram_f2.rsplit('/', 1)[-1])
        print('there are', cram_chromosomeHashTable_comparisons[2], 'reads in', cram_f2.rsplit('/', 1)[-1], 'not found mapped in', cram_f1.rsplit('/', 1)[-1])
        print('there are', cram_chromosomeHashTable_comparisons[3], 'reads in that were mapped differently between', cram_f2.rsplit('/', 1)[-1], ' and ', cram_f1.rsplit('/', 1)[-1])


final_time = time.time() - total_time
print('total exectuion time', round(final_time, 3), 'seconds (', round(final_time/60, 3), 'minutes).')
