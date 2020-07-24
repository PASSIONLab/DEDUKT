#!/usr/bin/python
#
# Requires python >= 3.6.4
#		change "...= match[match.lastindex]" for backwards compatibility...
#
#
# TODO: totalReads and totalBases are estimated not exact
#
#
from sys import argv
import re
import csv
import copy as cp

if (len(argv) < 2):
	print('At least 1 input file required for parsing. Exiting.')
	exit()

def search_fmt(pattern): return '('+pattern+')'

ints_pttn='[0-9]*'
ints_find=search_fmt(ints_pttn)
dec_pttn='[0-9]*\.*[0-9]*'
dec_find=search_fmt(dec_pttn)
file_pttn='[a-zA-Z0-9\_\-\.]*'
file_find=search_fmt(file_pttn)

ofileprefix=re.match('.*\-'+ints_find+'\..*', argv[1])
if (ofileprefix): ofileprefix=ofileprefix[ofileprefix.lastindex]
ofilesuffix=re.match('.*\-'+ints_find+'\..*', argv[len(argv)-1])
if (ofilesuffix): ofilesuffix=ofilesuffix[ofilesuffix.lastindex]
if (not ofileprefix and not ofilesuffix): ofileprefix=argv[1]
else: ofileprefix=ofileprefix+'-'+ofilesuffix

parsing_dict = {'Host Name: '+file_find:'Hostname',
				'Starting diBELLA version '+file_find+' on '+ints_pttn+' threads':'Version',
				'Running in mode: (.*)':'DbgOrRls_Mode',
				'Data: (.*)':'Data_shortname',
################Parallelism parameters
#				'Number of nodes: '+ints_find:'Num_nodes',
				'-n '+ints_find:'Num_total_threads',
				'srun .* -N '+ints_find+'|Number of nodes: '+ints_find:'Num_nodes',
				'-S '+ints_find:'Num_threads_per_numa',
				'-j '+ints_find:'Num_threads_per_cpu',
				'--cpus-per-proc '+ints_find:'Num_threads_per_cpu',
################Stages started/finished and optimizations run
				'# Starting stage loadfq .*':'Loadfq_started',
				'# Finished loadfq in .*':'Loadfq_finished',
				'# Starting stage kmermatch.*':'KmerMatch_started',
				'# Finished kmermatch.*':'KmerMatch_finished',
				'# Starting stage alignment.*':'Alignment_started',
				'# Finished alignment.*':'Alignment_finished',
				'Reading FASTQ file .* \(cached\)':'Bool_cachedIO',
				'Thread .* HeavyHitters performance: .*':'Bool_ranHH',
################Memory available/consumed
				'RANK 0: start_mem_free = '+dec_find:'MemGB_start_free_node_0',
				'# Memory remaining on node 0 after loadfq: '+dec_find+' GB':'MemGB_loadfq_end_node_0',
				'# Memory remaining on node 0 after kmermatch-'+ints_pttn+': '+dec_find+' GB':'MemGB_KM_end_node_0',
				'# Memory remaining on node 0 after alignment-'+ints_pttn+': '+dec_find+' GB':'MemGB_align_end_node_0',
################Runtime/quality parameters
				'-k '+ints_find:'Val_k',
				'Max k\-mer length \(internal\) = '+ints_find:'Val_max_internal_k',
				'-e '+ints_find:'Min_kmer_frequency',
				'-u '+ints_find:'Max_kmer_frequency',
				'-q '+ints_find:'Max_seeds_extended_per_pair',
				'Max alignments to attempt per overlapping read pair = '+ints_find:'Val_max_algns_per_overlap',
################loadfq Runtime statistics
				'# Starting stage loadfq -N '+ints_find+' .* ':'Num_Loadfq_threads_per_node',
				'Overall time for loadfq is '+dec_find+' s':'Time_Loadfq_overall',
################KmerMatch Runtime statistics
#				'kmermatch_main: .* totalReads: '+ints_find+' .*':'Total_reads', #TODO: this pattern retrieves estimated number of reads not exact number of reads...
#				'kmermatch_main: .* totalBases: '+ints_find:'Total_bases', #TODO: this pattern retrieves estimated number of reads not exact number of reads...
				'kmermatch_main: Pre 1st pass cardinality estimate, elapsed time: '+dec_find+' s':'Time_KM_cardinalityEst_elapsed',
				'ProcessFiles pass 1: Table size is: '+ints_pttn+' bits, '+dec_find+' MB':'SizeMB_BloomFilterLocal',
				'ProcessFiles pass 1: Optimal number of hash functions is : '+ints_find:'Num_BloomHashFunctions',
				'ProcessFiles pass 1: Average time taken for FASTQ reads is '+dec_find+'.*':'Time_KM_P1_FastqLoading_Avg',
				'ProcessFiles pass 1: Average time taken for packing reads is '+dec_find+'.*':'Time_KM_P1_Packing_Avg',
				'ProcessFiles pass 1: Average time taken for exchanging reads is '+dec_find+'.*':'Time_KM_P1_Exchanging_Avg',
				'ProcessFiles pass 1: Average time taken for processing reads is '+dec_find+'.*':'Time_KM_P1_LocalProcessing_Avg',
				'ProcessFiles pass 1: Average time taken for other FASTQ processing is '+dec_find+'.*':'Time_KM_P1_Other_Avg',
				'ProcessFiles pass 1: Average time taken for elapsed is '+dec_find+'.*':'Time_KM_P1_Elapsed_Avg',
				'ProcessFiles pass 2: Average time taken for FASTQ reads is '+dec_find+'.*':'Time_KM_P2_FastqLoading_Avg',
				'ProcessFiles pass 2: Average time taken for packing reads is '+dec_find+'.*':'Time_KM_P2_Packing_Avg',
				'ProcessFiles pass 2: Average time taken for exchanging reads is '+dec_find+'.*':'Time_KM_P2_Exchanging_Avg',
				'ProcessFiles pass 2: Average time taken for processing reads is '+dec_find+'.*':'Time_KM_P2_LocalProcessing_Avg',
				'ProcessFiles pass 2: Average time taken for other FASTQ processing is '+dec_find+'.*':'Time_KM_P2_Other_Avg',
				'ProcessFiles pass 2: Average time taken for elapsed is '+dec_find+'.*':'Time_KM_P2_Elapsed_Avg',
				'kmermatch_main: 2nd .* elapsed time: '+dec_find+' s':'Time_KM_P2_Elapsed',
				'kmermatch_main: Load imbalance for final k-mer counts: '+dec_find:'Val_KM_loadImbalance',
				'kmermatch_main: Total number of stored k-mers: '+ints_find:'Num_kmers_distinct_retained_global',
				'countTotalKmersAndCleanHash: Kmerscount hash includes '+ints_find+' distinct elements':'Num_kmers_distinct_afterBloom_global',
################Overlap Runtime statistics
				'extractAndExchangeOverlappingReads: Average time allocating per\-thread buffers before computing AxA\^T: '+dec_find+' s ':'Time_Overlap_AllocExch_Avg',
				'extractAndExchangeOverlappingReads: Average time for transforming local map \(kmer\-\>\[ReadId list\] to ReadId\-\>\[ReadId list\]\): '+dec_find+' s':'Time_Overlap_LocalAAT_Avg',
				'ReadOverlapper:exchangeReadPairs Average packing time: '+dec_find+' s':'Time_Overlap_IdPacking_Avg',
				'ReadOverlapper:exchangeReadPairs Average storage time: '+dec_find+' s':'Time_Overlap_IdStorage_Avg',
				'ReadOverlapper:exchangeReadPairs Min, max, avg and total \(sum\) time \(s\) for ReadId pair exchange: '+dec_find+' \(min\), '+dec_pttn+' \(max\), '+dec_pttn+' \(avg\), '+dec_pttn+' \(sum\)':'Time_Overlap_IdExch_Min',
				'ReadOverlapper:exchangeReadPairs Min, max, avg and total \(sum\) time \(s\) for ReadId pair exchange: '+dec_pttn+' \(min\), '+dec_find+' \(max\), '+dec_pttn+' \(avg\), '+dec_pttn+' \(sum\)':'Time_Overlap_IdExch_Max',
				'ReadOverlapper:exchangeReadPairs Min, max, avg and total \(sum\) time \(s\) for ReadId pair exchange: '+dec_pttn+' \(min\), '+dec_pttn+' \(max\), '+dec_find+' \(avg\), '+dec_pttn+' \(sum\)':'Time_Overlap_IdExch_Avg',
				'ReadOverlapper:exchangeReadPairs Min, max, avg and total \(sum\) time \(s\) for ReadId pair exchange: '+dec_pttn+' \(min\), '+dec_pttn+' \(max\), '+dec_pttn+' \(avg\), '+dec_find+' \(sum\)':'Time_Overlap_IdExch_Sum',
				'ReadOverlapper:exchangeReadNames : total \(sum\) and average \(packing\) times were '+dec_pttn+' s and '+dec_find+' s':'Time_Overlap_NamePacking_Avg',
				'ReadOverlapper:exchangeReadNames : total \(sum\) and average \(store\) times were '+dec_pttn+' s and '+dec_find+' s ':'Time_Overlap_NameStorage_Avg',
				'ReadOverlapper:exchangeReadNames : total \(sum\), average, min, and max \(exchange\) times \(s\) were: '+dec_find+', '+dec_pttn+', '+dec_pttn+', '+dec_pttn:'Time_Overlap_NameExch_Sum',
				'ReadOverlapper:exchangeReadNames : total \(sum\), average, min, and max \(exchange\) times \(s\) were: '+dec_pttn+', '+dec_find+', '+dec_pttn+', '+dec_pttn:'Time_Overlap_NameExch_Avg',
				'ReadOverlapper:exchangeReadNames : total \(sum\), average, min, and max \(exchange\) times \(s\) were: '+dec_pttn+', '+dec_pttn+', '+dec_find+', '+dec_pttn:'Time_Overlap_NameExch_Min',
				'ReadOverlapper:exchangeReadNames : total \(sum\), average, min, and max \(exchange\) times \(s\) were: '+dec_pttn+', '+dec_pttn+', '+dec_pttn+', '+dec_find:'Time_Overlap_NameExch_Max',
				'ReadOverlapper:outputReadOverlaps : max overlap\-map traversal time: '+dec_find+' s ':'Time_Overlap_LocalNameTraversal_Max',
				'ReadOverlapper:outputReadOverlaps : avg overlap\-map traversal time: '+dec_find+' s ':'Time_Overlap_LocalNameTraversal_Avg',
				'ReadOverlapper:outputReadOverlaps : max overlap output time: '+dec_find+' s ':'Time_Overlap_Output_Max',
				'ReadOverlapper:outputReadOverlaps : avg overlap output time: '+dec_find+' s ':'Time_Overlap_Output_Avg',
				'kmermatch_main: Total time computing, outputting, and cleaning\-up overlaps: '+dec_find+' s':'Time_Overlap_Elapsed',
				'kmermatch\-'+ints_pttn+': Overall time for kmermatch_main is '+dec_find+' s':'Time_KM_Elapsed',
################Alignment Data volumes
				'Total global number of overlaps: '+ints_find:'Num_Overlaps_Sum', #TODO
				'Aligner:computeLocalAlignments Local alignment maximum \(global elapsed\) time for '+ints_find:'Num_Algn_Attempted_Sum',
				'Aligner:computeLocalAlignments Min, max, and avg number of local alignments computed: '+ints_pttn+', '+ints_pttn+', '+ints_find:'Num_Algn_CompPP_Avg',
				'Aligner:computeLocalAlignments Min, max, and avg number of local alignments computed: '+ints_find:'Num_Algn_CompPP_Min',
				'Aligner:computeLocalAlignments Min, max, and avg number of local alignments computed: '+ints_pttn+', '+ints_find:'Num_Algn_CompPP_Max',
################Alignment Computation Time
				'Aligner:loadSequences Average time taken for FASTQ reads is '+dec_find:'Time_Algn_FastqLoad_Avg',
				'Aligner:loadSequences Average time taken for elapsed * is '+dec_find+', *':'Time_Algn_FastqLoad_Elapsed',
				'Aligner:computeLocalAlignments Local alignment maximum \(global elapsed\) time for '+ints_pttn+' total attempted alignments is '+dec_find:'Time_Algn_Comp_Global',
				'Aligner:computeLocalAlignments Min, max, avg and total \(sum\) time \(s\) for alignment computation: '+dec_pttn+', '+dec_pttn+', '+dec_find:'Time_Algn_Comp_Avg',
				'Aligner:computeLocalAlignments Min, max, avg and total \(sum\) time \(s\) for alignment computation: '+dec_find:'Time_Algn_Comp_Min',
				'Aligner:computeLocalAlignments Min, max, avg and total \(sum\) time \(s\) for alignment computation: '+dec_pttn+', '+dec_find:'Time_Algn_Comp_Max',
				'Aligner:computeLocalAlignments Min, max, avg and total \(sum\) time \(s\) for alignment computation: '+dec_pttn+', '+dec_pttn+', '+dec_pttn+', '+dec_find:'Time_Algn_Comp_Sum',
################Alignment Communication volumes
				'Aligner:exchangeSequences:iteration '+ints_pttn+': '+ints_find+' average reads requested.*':'Num_ReadExchReads_Avg',
				'Aligner:exchangeSequences:iteration '+ints_pttn+': '+ints_find+' min reads requested.*':'Num_ReadExchReads_Min',
				'Aligner:exchangeSequences:iteration '+ints_pttn+': '+ints_find+' max reads requested.*':'Num_ReadExchReads_Max',
				'Aligner:exchangeSequences:iteration '+ints_pttn+': .*?'+ints_find+' total reads requested':'Num_ReadExchReads_Sum',
				'Aligner:exchangeSequences:iteration '+ints_pttn+': sending '+ints_find+' min read bytes per rank':'Num_ReadExchSnd_bytes_Min',
				'Aligner:exchangeSequences:iteration '+ints_pttn+': sending '+ints_find+' max read bytes per rank':'Num_ReadExchSnd_bytes_Max',
				'Aligner:exchangeSequences:iteration '+ints_pttn+': receiving '+ints_find+' min read bytes per rank':'Num_ReadExchRcv_bytes_Min',
				'Aligner:exchangeSequences:iteration '+ints_pttn+': receiving '+ints_find+' max read bytes per rank':'Num_ReadExchRcv_bytes_Max',
################Alignment Communication Time
				'Aligner:exchangeSequences Sequence Exchange: total \(sum\) and average \(packing\) times were '+dec_find:'Time_Algn_ReadPacking_Sum',
				'Aligner:exchangeSequences Sequence Exchange: total \(sum\) and average \(packing\) times were '+dec_pttn+' s and '+dec_find:'Time_Algn_ReadPacking_Avg',
				'Aligner:exchangeSequences Sequence Exchange: total \(sum\) and average \(store\) times were '+dec_pttn+' s and '+dec_find:'Time_Algn_ReadStorage_Avg',
				'Aligner:exchangeSequences Sequence Exchange: total \(sum\) and average \(store\) times were '+dec_find: 'Time_Algn_ReadStorage_Sum',
				'Aligner:exchangeSequences Sequence Exchange: total \(sum\), average, min, and max \(exchange\) times \(s\) were: '+dec_pttn+', '+dec_find:'Time_ReadExch_Avg',
				'Aligner:exchangeSequences Sequence Exchange: total \(sum\), average, min, and max \(exchange\) times \(s\) were: '+dec_pttn+', '+dec_pttn+', '+dec_find:'Time_ReadExch_Min',
				'Aligner:exchangeSequences Sequence Exchange: total \(sum\), average, min, and max \(exchange\) times \(s\) were: '+dec_pttn+', '+dec_pttn+', '+dec_pttn+', '+dec_find:'Time_ReadExch_Max',
				'Aligner:exchangeSequences Sequence Exchange: total \(sum\), average, min, and max \(exchange\) times \(s\) were: '+dec_find:'Time_ReadExch_Sum',
################I/O Time
				'Loaded overlaps in '+dec_find+' s':'Time_Algn_OvlpLoading_Elapsed',
				'Aligner:computeLocalAlignments Min, max, avg and total \(sum\) time \(s\) for alignment IO: '+dec_pttn+', '+dec_pttn+', '+dec_find:'Time_Algn_IO_Avg',
				'Aligner:computeLocalAlignments Min, max, avg and total \(sum\) time \(s\) for alignment IO: '+dec_find:'Time_Algn_IO_Min',
				'Aligner:computeLocalAlignments Min, max, avg and total \(sum\) time \(s\) for alignment IO: '+dec_pttn+', '+dec_find:'Time_Algn_IO_Max',
				'Aligner:computeLocalAlignments Min, max, avg and total \(sum\) time \(s\) for alignment IO: '+dec_pttn+', '+dec_pttn+', '+dec_pttn+', '+dec_find:'Time_Algn_IO_Sum',
				'Aligner:computeLocalAlignments Local alignment maximum \(global elapsed\) time for '+ints_pttn+' total attempted alignments is '+dec_find:'Time_Algn_Elapsed',
################File and other handles
				'Input fastq file list: '+file_find:'Val_input_fastq_list_file',
				'Date: (.*$)':'Host_date',
				'Run directory is: (.*)':'Run_directory'
}


def writemycsv(fields, all_csv_rows, ofilename):
	with open(ofilename, 'w') as csvfile:
		writer = csv.DictWriter(csvfile, fieldnames=fields)
		writer.writeheader()
		for row in all_csv_rows:
			writer.writerow(row)


def get_par_pttn():
	return {'Num_nodes':'srun .* -N '+ints_find+'|[Nn]umber of [Nn]odes: '+ints_find,
			'Num_total_threads':'-n '+ints_find}

def get_exch_dict():
	return {'Num_nodes':'', 'Num_total_threads':'','ExchStep':'','TotalIters':'0','Iter':'',
		'MinSnd':'','MaxSnd':'', 'MinRcv':'', 'MaxRcv':'', 'MinTime':'', 'MaxTime':'', 'SourceFile':''}

def get_exch_headers(): return list(get_exch_dict().keys())

def get_hll_merge_dict():
	return {'Num_nodes':'', 'Num_total_threads':'', 'Data_shortname':'',
		'Time_Merge_Avg':'', 'Time_Merge_Max':'', 'SourceFile':''}

def get_hll_merge_headers(): return list(get_hll_merge_dict().keys())

def get_hll_dict():
	return {'Num_nodes':'', 'Num_total_threads':'', 'Data_shortname':'',
	'Val_FileName':'', 'Time_Fastq_Avg':'', 'Time_Fastq_Max':'', 'Time_InsertIntoHLL_Avg':'',
	'Time_InsertIntoHLL_Max':'', 'Time_Reduce_Avg':'', 'Time_Reduce_Max':'',
	'Time_Elapsed_Avg':'', 'Time_Elapsed_Max':'', 'SourceFile':'' }

def get_hll_headers(): return list(get_hll_dict().keys())

'''
Assumes parallelism features (e.g. number of nodes, ranks,...)
'''
def parse_results(texttosearch, sourcefile, dataname, pttns_dict, csv_dict, pttn_to_csv_mapping):
	par_pttns=get_par_pttn()
	pvals={p_label:re.search(p_pttn, texttosearch, re.MULTILINE)[1] for p_label,p_pttn in par_pttns.items()}
	matches={key:re.findall(pttn, texttosearch, re.MULTILINE) for (key,pttn) in pttns_dict.items()}
	csv_rows=[]
	for (key,vals) in matches.items():
		for v in vals:
			row=cp.deepcopy(csv_dict)
			for (plab,pval) in pvals.items(): row[plab]=pval
			row['SourceFile']=sourcefile
			row['Data_shortname']=dataname
			for (i,l) in pttn_to_csv_mapping: row[l]=v[i]
			csv_rows.append(row)
	return csv_rows

'''
Example:
[('human54x_p1.fastq', '0.158384', '0.207074', '57.2347', '57.882', '0.000941072', '1.44549', '57.4004', '58.0588')]
'''
def parse_hll_merge(texttosearch, sourcefile, dataname):
	merge_pttns={'hll_merge':'ProudlyParallelCardinalityEstimate:Total time taken for HH Merge '+dec_pttn+' avg HH merge: '+dec_find+' max: '+dec_find}
	csv_dict=get_hll_merge_dict()
	return parse_results(texttosearch, sourcefile, dataname, merge_pttns, csv_dict, list(zip([0,1], list(csv_dict.keys())[3:5])) )

'''
Example:
[('0.00475528', '0.0388298')]
'''
def parse_hll_summary(texttosearch, sourcefile, dataname):
	hll_pttns={'hll_summary':'ProudlyParallelCardinalityEstimate: Total time taken for FASTQ file '+file_find+' size: '+ints_pttn+' avg: '+dec_find+' max: '+dec_find+' min: '+dec_pttn+' HLL avg: '+dec_find+' HLL max: '+dec_find+' HH avg: '+dec_find+' HH max: '+dec_find+' elapsed_avg: '+dec_find+' elapsed_max: '+dec_find+' seconds'}
	csv_dict=get_hll_dict()
	return parse_results(texttosearch, sourcefile, dataname, hll_pttns, csv_dict, list(zip(list(range(0,9)), list(csv_dict.keys())[3:12])) )

# KmerMatch:Exchange exchange iteration 0 pass 1: avg: 0.938 MB at 0.938 MB/s, max: 0.939 MB at 0.939 MB/s, total: 0.704 GB in 0.939 s  [2,4,6] [3,5,7]
# ReadOverlapper:exchangeReadPairs exchange iteration 0, total 6080.028 MEGA ReadId pairs, Megabytes in ReadId(1) and Position(2) exchanges: avg1: 126.6672 MB, avg2: 31.6668 MB, max1: 126.9263 MB, max2: 31.7316 MB, total1: 95.0004 GB, total2: 23.7501 GB, and avg 34.3656 MB/s, max 32.9977 MB/s
# ReadOverlapper:exchangeReadNames exchange iteration 0, global total 1533.625 MEGA ReadIds exchanged, avg: 15.975 MB and 29.631 MB/s, max: 15.975 MB and 28.287 MB/s, total: 11.981 GB and 21.215 GB/s
# ReadOverlapper:exchangeReadNames exchange iteration 0, global avg: 155.710 MB and 67.5534 MB/s, max: 160.665 MB and 69.2845 MB/s, total: 116.782 GB and 50.3607 GB/s
def parse_exch_iters(texttosearch, sourcefile):
	par_pttns=get_par_pttn()
	exch_pttns={'kmermatch_p1_exch':'KmerMatch:Exchange exchange iteration '+ints_find+' pass 1: sent min '+ints_find+' bytes, sent max '+ints_find+' bytes, recv min '+ints_find+' bytes, recv max '+ints_find+' bytes, in min '+dec_find+' s, max '+dec_find+' s',
				'kmermatch_p2_exch':'KmerMatch:Exchange exchange iteration '+ints_find+' pass 2: sent min '+ints_find+' bytes, sent max '+ints_find+' bytes, recv min '+ints_find+' bytes, recv max '+ints_find+' bytes, in min '+dec_find+' s, max '+dec_find+' s',
				'overlap_rid_exch':'ReadOverlapper:exchangeReadPairs exchange iteration '+ints_find+', '+ints_find+' bytes per ReadId, '+ints_find+' bytes per Position: sent min '+ints_find+' IDs\/Positions, sent max '+ints_find+' IDs\/Positions, recv min '+ints_find+' IDs\/Positions, recv max '+ints_find+' IDs\/Positions, in min '+dec_find+' s, max '+dec_find+' s',
#				'labeling_rid_exch':'ReadOverlapper:exchangeReadNames exchange iteration '+ints_find+', global total '+dec_find+' MEGA ReadIds exchanged, avg: '+dec_find+' MB and '+dec_find+' MB/s, max: '+dec_find+' MB and '+dec_find+' MB/s, total: '+dec_find+' GB and '+dec_find+' GB/s',
#				'labeling_char_exch':'ReadOverlapper:exchangeReadNames exchange iteration '+ints_find+', global avg: '+dec_find+' MB and '+dec_find+' MB/s, max: '+dec_find+' MB and '+dec_find+' MB/s, total: '+dec_find+' GB and '+dec_find+' GB/s'
				}
	csv_dict=get_exch_dict()
	pvals={p_label:re.search(p_pttn, texttosearch, re.MULTILINE)[1] for p_label,p_pttn in par_pttns.items()}
	matches={key:re.findall(pttn, texttosearch, re.MULTILINE) for (key,pttn) in exch_pttns.items()}
	csv_rows=[]
	rid_size=0
	pos_size=0
	for (key,vals) in matches.items():
		for v in vals:
			row=cp.deepcopy(csv_dict)
			for (plab,pval) in pvals.items(): row[plab]=pval
			row['SourceFile']=sourcefile
			row['ExchStep']=key
			row['TotalIters']=len(vals)
			row['Iter']=v[0]
			keylist=[]
			if (key in ['kmermatch_p1_exch','kmermatch_p2_exch']):
				keylist=zip(range(1,7),list(csv_dict.keys())[5:11])
			if (key is 'overlap_rid_exch'):
				keylist=zip(range(3,10),list(csv_dict.keys())[5:11])
				rid_size=v[1]
				pos_size=v[2]
			for (i,l) in keylist: row[l]=v[i]
			csv_rows.append(row)
	print("sizeof(ReadId)="+str(rid_size)+", sizeof(PosInRead)="+str(pos_size))
	return csv_rows

all_csv_rows=[]
exch_csv_rows=[]
hll_csv_rows=[]
hll_merge_csv_rows=[]
for filename in argv[1:]:
	csv_row=dict()
	with open(filename,'rt') as infile:
		intext = infile.read()
# parse the per-iteration exchange information
		exchrows=parse_exch_iters(intext, filename)
		for row in exchrows: exch_csv_rows.append(row)
# parse the rest
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
	if(csv_row['Min_kmer_frequency'] == ''): csv_row['Min_kmer_frequency']='2' #default value with -e is not set
	if('Loadfq_started' in csv_row): csv_row['Loadfq_started']='True'
	else: csv_row['Loadfq_started']='False'
	if('Loadfq_finished' in csv_row): csv_row['Loadfq_finished']='True'
	else: csv_row['Loadfq_finished']='False'
	if('KmerMatch_started' in csv_row): csv_row['KmerMatch_started']='True'
	else: csv_row['KmerMatch_started']='False'
	if('KmerMatch_finished' in csv_row): csv_row['KmerMatch_finished']='True'
	else: csv_row['KmerMatch_finished']='False'
	if('Alignment_started' in csv_row): csv_row['Alignment_started']='True'
	else: csv_row['Alignment_started']='False'
	if('Alignment_finished' in csv_row): csv_row['Alignment_finished']='True'
	else: csv_row['Alignment_finished']='False'
	if('Bool_cachedIO' in csv_row): csv_row['Bool_cachedIO']='True'
	else: csv_row['Bool_cachedIO']='False'
	if('Bool_ranHH' in csv_row):
		csv_row['Bool_ranHH']='True'
#		if (csv_row['Data_shortname'] in csv_row):
		hllrows=parse_hll_summary(intext, filename, csv_row['Data_shortname'])
		for row in hllrows: hll_csv_rows.append(row)
		hllmerge_rows=parse_hll_merge(intext, filename, csv_row['Data_shortname'])
		for row in hllmerge_rows: hll_merge_csv_rows.append(row)
	else: csv_row['Bool_ranHH']='False'
	csv_row['Job_out_filename']=filename
	all_csv_rows.append(csv_row)

fields=[v for v in parsing_dict.values()]+['Job_out_filename']
writemycsv(fields, all_csv_rows, ofileprefix+'.csv')
writemycsv(get_exch_headers(), exch_csv_rows, ofileprefix+'_exchs.csv')
writemycsv(get_hll_headers(), hll_csv_rows, ofileprefix+'_hll.csv')
writemycsv(get_hll_merge_headers(), hll_merge_csv_rows, ofileprefix+'_hll_merge.csv')
