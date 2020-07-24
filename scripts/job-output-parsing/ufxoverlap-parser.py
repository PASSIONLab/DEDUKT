#!/usr/bin/python
#
# Requires python >= 3.6.4
#		change "...= match[match.lastindex]" for more backwards compatibility...
#
from sys import argv
import re
import csv

if (len(argv) < 2):
	print('One or more input files required for parsing.')
	exit()

'''
Sample input file:

Run directory is: diBELLA_64n_16c_ecoli30x_k17_e2_h7
Thu Mar  1 23:12:47 EST 2018
==================== Starting ufx ====================
aprun -n 1024 -N 16 -S 8 /lustre/atlas/scratch/mme/csc103/diBELLA-install/bin/ufx-8-32 -k 17 -i all_fastq.txt -e 2 -h 7
Thread 0 [INFO] Thu Mar  1 23:12:51 2018 [mach.h:273-(null)]: Found 8 cores_per_node
STAGE ufx_main -k 17 -i all_fastq.txt -e 2 -h 7 
HipMER k-mer characterization (ufx generation) step
You are running with the following settings:
K-mer length = 17
Max k-mer length (internal) = 32
Max read-records per k-mer = 7
pacbio_filtered.fastq : 266 MB
Estimated our cardinality: 139020602 totalReads: 15659 totalBases: 139256227
Stage 1 - Cardinality: 0.067 s
Table size is: 761857 bits, 0.0908204 MB
Optimal number of hash functions is : 5
Reading FASTQ file pacbio_filtered.fastq 
Cardinality for pacbio_filtered.fastq: 139020602
active ranks moreSeqs: 1024 moreToExchange: 0 moreFiles: 0, rank 0 moreSeqs: 1 moreToExchange: 0 moreFiles: 0 pack_time: 0.191 exchange_time: 1.936 process_time: 0.000 elapsed: 11.848
Average time taken for FASTQ reads is 0.002, myelapsed 0.004
Average time taken for packing reads is 0.201, myelapsed 0.192
Average time taken for exchanging reads is 3.703, myelapsed 1.951
Average time taken for processing reads is 0.312, myelapsed 0.000
Average time taken for elapsed is 11.868, myelapsed 11.880
PASS 1
Read/distributed/processed reads of  ALL files  in 11.880 seconds
Stage 2 - Exchange and Insert: 11.881
Reading FASTQ file pacbio_filtered.fastq 
active ranks moreSeqs: 1024 moreToExchange: 0 moreFiles: 0, rank 0 moreSeqs: 1 moreToExchange: 0 moreFiles: 0 pack_time: 0.264 exchange_time: 1.918 process_time: 0.000 elapsed: 6.641
Average time taken for FASTQ reads is 0.002, myelapsed 0.004
Average time taken for packing reads is 0.257, myelapsed 0.264
Average time taken for exchanging reads is 1.972, myelapsed 1.922
Average time taken for processing reads is 0.390, myelapsed 0.000
Average time taken for elapsed is 6.653, myelapsed 6.659
PASS 2
Read/distributed/processed reads of  ALL files  in 6.659 seconds
Counting finished 
Kmerscount hash includes 28579790 distinct elements
Kmerscount non error kmers count is 51248739
Global max count is 27186
Large count histogram is of size 300000
Thread 0 [INFO] Thu Mar  1 23:13:10 2018 [UFXextended.cpp:1130-(null)]: Skipping heavy hitter counts
Erroneous count < 2 cases removed 
Kmerscount hash includes 10003311 distinct elements
Kmerscount non error kmers count is 29250345
Generated histograms
Stage 3 Count 7.827 s
Total count is 10003311
Load imbalance for final k-mer counts is 9.539
CardinalityEstimate 4.023 MEGA k-mers per sec/proc in 0.067 seconds
Bloom filter + hash pre-allocation  0.023 MEGA k-mers per sec/proc in 11.881 seconds
Actual hash-based counting  0.035 MEGA k-mers per sec/proc in 7.827 seconds
Stage 4 - loadImbalance: 0.050 s
Writing to binary via MPI-IO
Thread 0 [INFO] Thu Mar  1 23:13:11 2018 [UFXextended.cpp:1558-(null)]: Writing output in 1 batches of about 9768 per rank
Stage 6 - FileIO: 1.869 s
File write completed
Finished writing, here is the top of processor 0's data
Stage 7 - shmem storage: 0.000 s
Overall time for ufx-8-32 is 21.69 s
Application 16913124 resources: utime ~21568s, stime ~558s, Rss ~235516, inblocks ~26107382, outblocks ~4303955
==================== Starting overlap ====================
aprun -n 1024 -N 16 -S 8 /lustre/atlas/scratch/mme/csc103/diBELLA-install/bin/overlap-8-32 -k 17 -i all_fastq.txt.ufx.bin -h 7
Thread 0 [INFO] Thu Mar  1 23:13:22 2018 [mach.h:273-(null)]: Found 8 cores_per_node
starting overlap_main
parsing k 
k is 17
parsing i 
finished parsing i 
diBELLA read-to-read overlap step
You are running with the following settings:
Input kmer-to-reads file = all_fastq.txt.ufx.bin
K-mer length (minimum contiguous overlap for any read pairing) = 17
Max k-mer length (internal) = 32
Reliable k-mer max = 7
cached_io is 0
successfully initialized endReadId[] 
kmer_length 17
Filesize is 960317856 bytes, and number of records is 10003311
Thread 0: Allocating memory for 9768 reads: 312576 bytes
Threads done with I/O
Average ufx loading time: 0.173 s 
Average hash map insertion (building) time: 0.043 s 
Average time for transforming local map (kmer->[ReadId list] to ReadId->[ReadId list]): 0.06 s 
Average packing time: 0.01 s 
Average storage time: 0.07 s 
Min, max, avg and total (sum) time (s) for ReadId pair exchange: 9.600 (min), 9.637 (max), 9.624 (avg), 9854.618191 (sum)
Additional I/O time: 0.264 s 
Generated read-overlap histograms
Application 16913131 resources: utime ~10723s, stime ~521s, Rss ~24768, inblocks ~13523789, outblocks ~3108241
Thu Mar  1 23:13:38 EST 2018

Selected fields:

UFX
Stage 2 - Exchange and Insert: 11.881
Stage 3 Count 7.827 s
Load imbalance for final k-mer counts is 9.539
Stage 6 - FileIO: 1.869 s
Overall time for ufx-8-32 is 21.69 s

READ OVERLAP
Average time for transforming local map (kmer->[ReadId list] to ReadId->[ReadId list]): 0.06 s 
Average packing time: 0.01 s 
Average storage time: 0.07 s 
Min, max, avg and total (sum) time (s) for ReadId pair exchange: 9.600 (min), 9.637 (max), 9.624 (avg), 9854.618191 (sum)
Additional I/O time: 0.264 s 
'''

ofilename=argv[1]+'.csv'

def search_fmt(pattern): return '('+pattern+')'

ints_pttn='[0-9]*'
ints_find=search_fmt(ints_pttn)
dec_pttn='[0-9]*\.*[0-9]*'
dec_find=search_fmt(dec_pttn)
file_pttn='[a-zA-Z0-9_\-\.]*'
file_find=search_fmt(file_pttn)

parsing_dict = {'Host Name: '+file_find:'Hostname', # TODO start echo-ing $HOSTNAME in job script
################Runtime parameters
				'Data: (.*)':'Data_shortname',
				'Running in mode: (.*)':'DbgOrRls_Mode',
				'Number of nodes: '+ints_find:'Num_nodes',
				'-n '+ints_find:'Num_total_mpi_tasks',
				'-N '+ints_find:'Num_mpi_tasks_per_node',
				'-S '+ints_find:'Num_mpi_tasks_per_numa',
				'-j '+ints_find:'Num_tasks_per_processor',
				'-k '+ints_find:'Val_k',
				'Max k\-mer length \(internal\) = '+ints_find:'Val_internal_k_max',
				'-e '+ints_find:'Val_reliable_min',
				'-h '+ints_find:'Val_reliable_max',
################UFX Runtime statistics
				'Stage 2 \- Exchange and Insert: '+dec_find:'Time_ufx_stage2_exchange_insert',
				'Stage 3 Count '+dec_find+' s':'Time_ufx_stage3_count',
				'Load imbalance for final k\-mer counts is '+dec_find:'Num_ufx_load_imbalance_kmercounts',
				'Stage 6 \- FileIO: '+dec_find+' s':'Time_ufx_stage6_io',
				'Overall time for '+file_pttn+' is '+dec_find+' s':'Time_ufx_overall',
################Overlap Runtime statistics				
				'Average time for transforming local map \(kmer\-\>\[ReadId list\] to ReadId\-\>\[ReadId list\]\): '+dec_find+' s':'Time_avg_local_AAT',
				'Average packing time: '+dec_find+' s':'Time_avg_overlap_packing',
				'Average storage time: '+dec_find+' s':'Time_avg_overlap_storage',
				'Min, max, avg and total \(sum\) time \(s\) for ReadId pair exchange: '+dec_find+' \(min\), '+dec_pttn+' \(max\), '+dec_pttn+' \(avg\), '+dec_pttn+' \(sum\)':'Time_min_exchange_overlaps',
				'Min, max, avg and total \(sum\) time \(s\) for ReadId pair exchange: '+dec_pttn+' \(min\), '+dec_find+' \(max\), '+dec_pttn+' \(avg\), '+dec_pttn+' \(sum\)':'Time_max_exchange_overlaps',
				'Min, max, avg and total \(sum\) time \(s\) for ReadId pair exchange: '+dec_pttn+' \(min\), '+dec_pttn+' \(max\), '+dec_find+' \(avg\), '+dec_pttn+' \(sum\)':'Time_avg_exchange_overlaps',
				'Min, max, avg and total \(sum\) time \(s\) for ReadId pair exchange: '+dec_pttn+' \(min\), '+dec_pttn+' \(max\), '+dec_pttn+' \(avg\), '+dec_find+' \(sum\)':'Time_sum_exchange_overlaps',
				'Additional I/O time: '+dec_find+' s':'Time_overlap_io',				 
################File and other handles				
				'Input kmer\-to\-reads file = '+file_find+'\.ufx\.bin':'Val_input_fastq_list_file',
				'Run directory is: (.*)':'Run_directory',
				'Date: (.*$)':'Host_date'
}

all_csv_rows=[]
for filename in argv[1:]:
	csv_row=dict()
	with open(filename,'rt') as infile:
		intext = infile.read()
		for (pattern,field_name) in parsing_dict.items():
#			print ('pattern: '+pattern)
			match = re.search(pattern, intext, re.MULTILINE)
			if(match):
				if(match.lastindex):
					csv_row[ parsing_dict[pattern] ] = match[match.lastindex]
#					print ('match[match.lastindex]: '+match[match.lastindex])
				else:
					csv_row[ parsing_dict[pattern] ] = ''
#					print('no such group')
	if(csv_row['Val_reliable_min'] == ''): csv_row['Val_reliable_min']='2' #default value with -e is not set
	csv_row['Job_out_filename']=filename
	all_csv_rows.append(csv_row)

with open(ofilename, 'w') as csvfile:
	fields=[v for v in parsing_dict.values()]+['Job_out_filename']
	writer = csv.DictWriter(csvfile, fieldnames=fields)
	writer.writeheader()
	for row in all_csv_rows:
		writer.writerow(row)
