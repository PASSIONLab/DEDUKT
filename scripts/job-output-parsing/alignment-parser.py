#!/usr/bin/python
#
# Requires python >= 3.6.4
#		change "...= match[match.lastindex]" for more backwards compatibility...
#
from sys import argv
import re
import csv

if (len(argv) < 2):
	print('At least 1 input file required for parsing. Exiting.')
	exit()

'''
Sample input file:

Thu Feb 22 16:52:58 EST 2018
==================== Starting ufx ====================
aprun -n 512 -N 8 -S 4 -j 1 /lustre/atlas/scratch/mme/csc103/diBELLA-install/bin/ufx-8-32 -k 17 -i all_fastq.txt -e 2 -h 7
==================== Starting overlap ====================
aprun -n 512 -N 8 -S 4 -j 1 /lustre/atlas/scratch/mme/csc103/diBELLA-install/bin/overlap-8-32 -k 17 -i all_fastq.txt.ufx.bin -h 7
==================== Starting alignment ====================
aprun -n 512 -N 8 -S 4 -j 1 /lustre/atlas/scratch/mme/csc103/diBELLA-install/bin/alignment-8-32 -k 17 -i all_fastq.txt -m overlaps-17
Thread 0 [INFO] Thu Feb 22 16:53:01 2018 [mach.h:273-(null)]: Found 8 cores_per_node
diBELLA read alignment step
You are running with the following settings:
K-mer (minimum seed) length = 17
Max k-mer length (internal) = 32
Max alignments to attempt per overlapping read pair = 65535
input fastq file list: all_fastq.txt
input overlap file: overlaps-17
output file name: alignments-17
pacbio_filtered.fastq : 266 MB
Reading FASTQ file pacbio_filtered.fastq 
Average time taken for FASTQ reads is 0.00399497, myelapsed 0.00612497
Average time taken for elapsed     is 0, myelapsed 1.87219
Loaded reads of  ALL files  in 1.87219 seconds
Loaded overlaps in 1.28 s
Finished sequence exchange 
Sequence Exchange: total (sum) and average (packing) times were 49.09 s and 0.10 s 
Sequence Exchange: total (sum) and average (store) times were 89.96 s and 0.18 s 
Sequence Exchange: total (sum), average, min, and max (exchange) times (s) were: 668.63, 1.31, 1.27, 1.35 
Local alignment maximum (global elapsed) time for 122683600 total attempted alignments is 2353.678
Min, max, and avg number of local alignments computed: 236352, 243253, 239616 
Min, max, avg and total (sum) time (s) for alignment computation: 1955.730, 2351.849, 2159.091, 1105454.440 
Min, max, avg and total (sum) time (s) for alignment IO: 1.606, 3.185, 2.324, 1189.857 
Application 16864545 resources: utime ~1205344s, stime ~919s, Rss ~397288, inblocks ~10285736, outblocks ~2969504
Thu Feb 22 17:32:22 EST 2018
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
				'K\-mer \(minimum seed\) length = '+ints_find:'Val_k',
				'Max k\-mer length \(internal\) = '+ints_find:'Val_internal_k_max',
				'-e '+ints_find:'Val_reliable_min',
				'-h '+ints_find:'Val_reliable_max',
				'Max alignments to attempt per overlapping read pair = '+ints_find:'Val_max_algn_limit_per_overlap',
################Data volumes
				'Total global number of overlaps: '+ints_find:'Num_sum_overlaps',
				'Local alignment maximum \(global elapsed\) time for '+ints_find:'Num_sum_attempted_alignments',
				'Min, max, and avg number of local alignments computed: '+ints_pttn+', '+ints_pttn+', '+ints_find:'Num_avg_algn_computed_pp',
				'Min, max, and avg number of local alignments computed: '+ints_find:'Num_min_algn_computed_pp',
				'Min, max, and avg number of local alignments computed: '+ints_pttn+', '+ints_find:'Num_max_algn_computed_pp',
################Computation Time				
				'Average time taken for FASTQ reads is '+dec_find:'Time_avg_fastq_loading',
				'Average time taken for elapsed * is '+dec_find:'Time_elapsed_fastq_loading',
				'Loaded overlaps in '+dec_find+' s':'Time_elapsed_overlap_loading',
				'Local alignment maximum \(global elapsed\) time for '+ints_pttn+' total attempted alignments is '+dec_find:'Time_global_algn_computation',
				'Min, max, avg and total \(sum\) time \(s\) for alignment computation: '+dec_pttn+', '+dec_pttn+', '+dec_find:'Time_avg_algn_computation',
				'Min, max, avg and total \(sum\) time \(s\) for alignment computation: '+dec_find:'Time_min_algn_computation',
				'Min, max, avg and total \(sum\) time \(s\) for alignment computation: '+dec_pttn+', '+dec_find:'Time_max_algn_computation',
				'Min, max, avg and total \(sum\) time \(s\) for alignment computation: '+dec_pttn+', '+dec_pttn+', '+dec_pttn+', '+dec_find:'Time_sum_algn_computation',
################Communication Time				
				'Sequence Exchange: total \(sum\) and average \(packing\) times were '+dec_find:'Time_sum_read_packing',
				'Sequence Exchange: total \(sum\) and average \(packing\) times were '+dec_pttn+' s and '+dec_find:'Time_avg_read_packing',
				'Sequence Exchange: total \(sum\) and average \(store\) times were '+dec_pttn+' s and '+dec_find:'Time_avg_read_storage',
				'Sequence Exchange: total \(sum\) and average \(store\) times were '+dec_find: 'Time_sum_read_storage',
				'Sequence Exchange: total \(sum\), average, min, and max \(exchange\) times \(s\) were: '+dec_pttn+', '+dec_find:'Time_avg_read_exchange',
				'Sequence Exchange: total \(sum\), average, min, and max \(exchange\) times \(s\) were: '+dec_pttn+', '+dec_pttn+', '+dec_find:'Time_min_read_exchange',
				'Sequence Exchange: total \(sum\), average, min, and max \(exchange\) times \(s\) were: '+dec_pttn+', '+dec_pttn+', '+dec_pttn+', '+dec_find:'Time_max_read_exchange',
				'Sequence Exchange: total \(sum\), average, min, and max \(exchange\) times \(s\) were: '+dec_find:'Time_sum_read_exchange',
################I/O Time	
				'Min, max, avg and total \(sum\) time \(s\) for alignment IO: '+dec_pttn+', '+dec_pttn+', '+dec_find:'Time_avg_algn_io',
				'Min, max, avg and total \(sum\) time \(s\) for alignment IO: '+dec_find:'Time_min_algn_io',
				'Min, max, avg and total \(sum\) time \(s\) for alignment IO: '+dec_pttn+', '+dec_find:'Time_max_algn_io',
				'Min, max, avg and total \(sum\) time \(s\) for alignment IO: '+dec_pttn+', '+dec_pttn+', '+dec_pttn+', '+dec_find:'Time_sum_algn_io',
################File and other handles				
				'input fastq file list: '+file_find:'Val_input_fastq_list_file',
				'input overlap file: '+file_find:'Val_input_overlap_file',
				'output file name: '+file_find:'Val_output_alignments_file',
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
