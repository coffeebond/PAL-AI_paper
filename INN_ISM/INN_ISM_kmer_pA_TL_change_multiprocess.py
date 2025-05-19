# version 20250321
import sys, subprocess, numpy, time, math, concurrent.futures, argparse, itertools, h5py, re, os
import numpy as np
from time import time
from datetime import datetime

# parameters for the input file
kmer_len_min = 3

chunk = 1000000 # every number of lines to prompt the progress
chunk_lines = 50000 # number of lines to be process for each core
n_threads = 20 # number of cores to use

###------functions------------------------------------------------------------------------
def timer(): # calculate runtime
	temp = str(time()-t_start).split('.')[0]
	temp =  '\t' + temp + 's passed...' + '\t' + str(datetime.now())
	return temp

def generate_kmers(kmer_len_min, kmer_len_max):
    #Generates all kmers
    kmers = []
    for l in range(kmer_len_min, kmer_len_max+1):
        kmers += ["".join(kmer) for kmer in list(itertools.product(["A","C","G","T"], repeat=l))]
    return kmers

def pwrite(f, text, timestamp=None):
	"""Write printed output to a file and print it to console."""
	log_text = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - {text}" if timestamp else text
	if f:
		f.write(log_text + '\n')
		f.flush()
	print(log_text)

def cn_convert(cn):
	if cn in range(4):
		if cn == 0:
			return('A')
		elif cn == 1:
			return('C')
		elif cn == 2:
			return('G')
		else:
			return('T')
	elif cn in list('ACGT'):
		if cn == 'A':
			return(0)
		elif cn == 'C':
			return(1)
		elif cn == 'G':
			return(2)
		else:
			return(3)
	else:
		print("Warning! Unrecogenized input: " + str(cn))
		return(None)

def exam_kmers_mp(idx_lst):
	gain_dtemp = {}
	loss_dtemp = {}
	global seq_ary, tl_ary, utr_len_limit
	for idx in idx_lst:
		if type(seq_ary[idx]) == bytes:
			seq = seq_ary[idx].decode('utf-8').upper()
		else:
			seq = seq_ary[idx].upper()
		tl_diff = tl_ary[idx]

		for pos in range(len(seq)):
			# only examine UTR positions
			if ('utr_bool' in globals() and utr_bool[pos] == 1) or ('utr_bool' not in globals()):

				# must be "ACGT" and withine UTR length limit to examine
				if seq[pos] in 'ACGT' and (pos + utr_len_limit) >= len(seq):
					# get all kmers covering the nucleotide of interest at this position
					loss_kmers = np.unique([seq[i:(i+l)] for i in range(max(pos-args.l+1,0), min(pos+1, len(seq)),1) for l in range(kmer_len_min, args.l+1, 1) if (i+l)>pos and (i+l)<=len(seq)])
					tl_loss = np.sum(tl_diff[pos,:])/3 # average of all mutations for calculating differnce due to loss of the kmer
					
					for kmer in loss_kmers:
						if re.search(r'[^ACGT]', kmer) is None: # make sure it only contains ACGT
							# update loss of kmer with running average and distance 
							if kmer in loss_dtemp:
								loss_dtemp[kmer][0] += 1
								avg_last = loss_dtemp[kmer][1]
								loss_dtemp[kmer][1] = avg_last + (tl_loss - avg_last) / loss_dtemp[kmer][0]
								loss_dtemp[kmer][2] = loss_dtemp[kmer][2] + (tl_loss - loss_dtemp[kmer][1]) * (tl_loss - avg_last)
							else:
								loss_dtemp[kmer] = np.asarray([1, tl_loss, 0])

					for m in range(4):
						if m != cn_convert(seq[pos]): # each of the three possible mutations
							tl_gain = tl_diff[pos, m]

							# get all kmers covering the mutated nucleotide of interest at this position
							new_seq = seq[:pos] + cn_convert(m) + seq[(pos+1):]
							gain_kmers = np.unique([new_seq[i:(i+l)] for i in range(max(pos-args.l+1,0), min(pos+1, len(seq)),1) for l in range(kmer_len_min, args.l+1, 1) if (i+l)>pos and (i+l)<=len(seq)])
							
							# remove kmers that have 'N'
							gain_kmers = [kmer for kmer in gain_kmers if re.search(r'[^ACGT]', kmer) is None]

							# update gain of kmers with running average and distance
							for kmer in gain_kmers:
								if kmer in gain_dtemp:
									gain_dtemp[kmer][0] += 1
									avg_last = gain_dtemp[kmer][1]
									gain_dtemp[kmer][1] = avg_last + (tl_gain - avg_last) / gain_dtemp[kmer][0]
									gain_dtemp[kmer][2] = gain_dtemp[kmer][2] + (tl_gain - gain_dtemp[kmer][1]) * (tl_gain - avg_last)
								else:
									gain_dtemp[kmer] = np.asarray([1, tl_gain, 0])
	return(gain_dtemp, loss_dtemp)
		
def add_to_first_dict(d1, d2):
	# d1 must have all keys that d2 have
	# https://math.stackexchange.com/questions/2971315/how-do-i-combine-standard-deviations-of-two-groups
	for key in d2:
		if key in d1:
			count_last = d1[key][0]
			d1[key][0] = d1[key][0] + d2[key][0]
			avg_last = d1[key][1]
			d1[key][1] = (count_last * d1[key][1] + d2[key][0] * d2[key][1]) / (count_last + d2[key][0])
			d1[key][2] = d1[key][2] + d2[key][2] + count_last * d2[key][0] * (avg_last - d2[key][1]) * (avg_last - d2[key][1]) / (count_last + d2[key][0])
		else: 
			sys.exit('No key ' + key + ' in dictionary when excuting function "add_to_first_dict"!\n')

def remove_pseudo_count(d):
	for k in d:
		if d[k][0] > 1:
			count_last = d[k][0]
			d[k][0] = d[k][0] - 1
			avg_last = d[k][1]
			d[k][1] = (d[k][1] * count_last - 0) / d[k][0]
			d[k][2] = d[k][2] - 0 - count_last / d[k][0] * avg_last * avg_last		
		else:
			d[k] = [0, 0, 0]

#####################################################################################################################
###------the script runs from there-------------------------------------------------------
t_start = time() # timer start

# parse the input
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest = 'i', type = str, help = 'input file in h5 format with CNN ISM result')
parser.add_argument('-r', '--reference', dest = 'r', type = str, help = 'reference poly(A) site file')
parser.add_argument('-l', '--length', dest = 'l', type = int, default = 8, help = 'longest kmer to examine')
parser.add_argument('-u', '--utr_length', dest = 'u', type = int, help = 'longest utr region to examine (from 3\' ends)')
parser.add_argument('-o', '--output', dest = 'o', type = str, help = 'output folder')
args = parser.parse_args()  

# output folder
if not os.path.exists(args.o):
	os.makedirs(args.o)
if not args.o.endswith('/'):
	args.o = args.o + '/'

idf = args.o + args.i.split('/')[-1].split('.')[0]
if args.u:
	idf = idf + f"_utr_len_limit_{args.u}"

f_log = open(idf + '_kmer_pA_weight_avg_dis_log.txt', 'w')
pwrite(f_log, 'Input file: ' + args.i)
pwrite(f_log, 'Output data in this folder:\n' + args.o)
pwrite(f_log, '\nExamining kmers with length from ' + str(kmer_len_min) + ' to ' + str(args.l) + '...' + timer())

# make a dictionary of poly(A) sites for filtering those are uncertain (low confidence)
pa_dict = {}
with open(args.r, 'r') as f:
	cols = f.readline().strip().split('\t')
	n_col_pa_id = cols.index('pA_id')
	n_col_confidence = cols.index('confidence')
	for line in f:
		lst = line.strip().split('\t')
		if lst[n_col_confidence] == 'high':
			if lst[n_col_pa_id] not in pa_dict:
				pa_dict[lst[n_col_pa_id]] = True
			else:
				sys.exit('Duplicated pA id found: ' + lst[n_col_pa_id])
pwrite(f_log, 'Number of confident pA sites in the reference file: ' + str(len(pa_dict)))

# make two master dictionary with kmer as keys, one is for gaining of a kmer and the other is for losing a kmer
# {kmer: [count of sequences, avg TL, dis TL}
# "avg" is the mean and "dis" is std^2 * (n-1)
# it's necessary to initiate with pseudo count at all position, otherwise dictionary combining in "add_to_first_dict" will have a "divide by zero error"
dict_gain = {} 
dict_loss = {}
for kmer in generate_kmers(kmer_len_min = kmer_len_min, kmer_len_max = args.l):
	dict_gain[kmer] = np.zeros(3)
	dict_gain[kmer][0] = dict_gain[kmer][0] + 1
	dict_loss[kmer] = np.zeros(3)
	dict_loss[kmer][0] = dict_loss[kmer][0] + 1

with h5py.File(args.i, 'r') as f:
	# get the poly(A) site ID
	if type(f['id'][:][0]) == bytes:
		idx_sele_lst = [i for i in range(f['id'][:].shape[0]) if f['id'][:][i].decode('utf-8') in pa_dict]
	else:
		idx_sele_lst = [i for i in range(f['id'][:].shape[0]) if f['id'][:][i] in pa_dict]

	pwrite(f_log, '\nNumber of confident pA sites in the input file: ' + str(len(idx_sele_lst)))

	# get the sequences
	seq_ary = f['seq'][:]

	# get the utr boolean label if it exists
	if 'utr_bool' in f:
		utr_bool = f['utr_bool'][:]

	# get the scores
	if 'diff_scores' in f.keys():
		tl_ary = f['diff_scores'][:]
	else:
		tl_ary = f['all_scores'][:] - f['wt_scores'][:][:, None, None]

	# wild-type nucleotides have NaN in the array. Replace NaN with 0.
	tl_ary[np.isnan(tl_ary)] = 0

	# set the maximal length of 3' UTR region to examine (from 3' ends)
	if args.u:
		utr_len_limit = max(args.u, 100)
	else:
		utr_len_limit = tl_ary.shape[1]
	pwrite(f_log, f"\nLength of the 3\' UTR region to examine (from the 3\' ends): {utr_len_limit}")

	with concurrent.futures.ProcessPoolExecutor(n_threads) as pool:
		futures = pool.map(exam_kmers_mp, np.array_split(range(len(idx_sele_lst)), n_threads))
		for d_gain, d_loss in futures: # combine data from outputs from all processes
			add_to_first_dict(dict_gain, d_gain)
			add_to_first_dict(dict_loss, d_loss)
	pwrite(f_log, 'Finished examining kmers...')				


pwrite(f_log, '\nRemove pseudo counts and output data...' + timer())
# need to remove the original pseudo count, which all have average weight 0 and standard deviation 0
remove_pseudo_count(dict_gain)
remove_pseudo_count(dict_loss)

with open(idf + '_kmer_gain_pA_weight_avg_dis.txt', 'w') as f:
	for kmer in sorted(dict_gain.keys(), key = lambda x: (len(x), x)):
		f.write(kmer + '\t' + '\t'.join(list(map(str, dict_gain[kmer]))) + '\n')

with open(idf + '_kmer_loss_pA_weight_avg_dis.txt', 'w') as f:
	for kmer in sorted(dict_loss.keys(), key = lambda x: (len(x), x)):
		f.write(kmer + '\t' + '\t'.join(list(map(str, dict_loss[kmer]))) + '\n')

pwrite(f_log, 'Final: ' + timer())		
f_log.close()
 
