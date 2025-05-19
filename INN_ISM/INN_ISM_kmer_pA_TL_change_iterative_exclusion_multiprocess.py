# version 20250330
import sys, subprocess, time, math, concurrent.futures, argparse, itertools, h5py, re, os
import numpy as np
import pandas as pd
from time import time
from datetime import datetime
from scipy import stats

chunk = 1000000 # every number of lines to prompt the progress
chunk_lines = 50000 # number of lines to be process for each core
n_threads = 40 # number of cores to use

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

def generate_mut_seqs_by_pos(seq, pos):
	if pos not in range(len(seq)):
		return([])

	original_nt = seq[pos]
	new_seq_lst = []
	for new_nt in [nt for nt in ['A', 'C', 'G', 'T'] if nt != original_nt]:
		new_seq_lst.append(seq[:pos] + new_nt + seq[pos+1:])
	return(new_seq_lst)

def exam_kmers_mp(idx_lst):
	gain_dtemp = {}
	loss_dtemp = {}
	global seq_ary, tl_ary, utr_len_limit, kmer_excl_loss, kmer_excl_gain, args
	
	for idx in idx_lst:
		if type(seq_ary[idx]) == bytes:
			seq = seq_ary[idx].decode('utf-8').upper()
		else:
			seq = seq_ary[idx].upper()
		tl_diff = tl_ary[idx] 

		for pos in range(len(seq)):

			# only examine UTR positions
			if ('utr_bool' in globals() and utr_bool[pos] == 1) or ('utr_bool' not in globals()):

				# must be "ACGT" and within UTR length limit to examine
				if seq[pos] in 'ACGT' and (pos + utr_len_limit) >= len(seq):
					
					#--- this part is for loss of kmers
					
					# get all kmers covering the nucleotide of interest at this position
					loss_kmers = set()
					for i in range(pos-args.l+1, pos+1): # iterate through all relevant positions
						if (i >= 0) and (i + args.l <= len(seq)): 
							original_kmer = seq[i:(i+args.l)]

							# all loss kmers covering this
							if re.search(r'[^ACGT]', original_kmer) is None: # make sure it only contains ACGT	
								loss_kmers.add(original_kmer)

					if len(loss_kmers) != 0 and (not loss_kmers & kmer_excl_loss): # none of kmers lost due to the mutation at this position is in the exclusion list
						for i in range(pos-args.l+1, pos+1): # iterate through all relevant positions
							if (i >= 0) and (i + args.l <= len(seq)): 
								original_kmer = seq[i:(i+args.l)]

								if re.search(r'[^ACGT]', original_kmer) is None: # make sure it only contains ACGT
									nt_pos_in_kmer = pos - i
									
									# update loss of original_kmer with running average and distance 
									nt_idx = np.array([idx for idx in range(4) if idx != cn_convert(seq[pos])])
									if original_kmer in loss_dtemp:
										loss_dtemp[original_kmer][nt_idx, nt_pos_in_kmer, 0] += 1
										avg_last = np.array(loss_dtemp[original_kmer][nt_idx, nt_pos_in_kmer, 1])
										loss_dtemp[original_kmer][nt_idx, nt_pos_in_kmer, 1] = avg_last + (tl_diff[pos,nt_idx] - avg_last) / loss_dtemp[original_kmer][nt_idx, nt_pos_in_kmer, 0]
										loss_dtemp[original_kmer][nt_idx, nt_pos_in_kmer, 2] = loss_dtemp[original_kmer][nt_idx, nt_pos_in_kmer, 2] + (tl_diff[pos,nt_idx] - loss_dtemp[original_kmer][nt_idx, nt_pos_in_kmer, 1]) * (tl_diff[pos,nt_idx] - avg_last)
									else:
										loss_dtemp[original_kmer] = np.zeros((4, len(original_kmer), 3))
										loss_dtemp[original_kmer][~reverse_one_hot(original_kmer),:] = np.nan
										loss_dtemp[original_kmer][nt_idx, nt_pos_in_kmer, 0] = 1
										loss_dtemp[original_kmer][nt_idx, nt_pos_in_kmer, 1] = tl_diff[pos,nt_idx]
										loss_dtemp[original_kmer][nt_idx, nt_pos_in_kmer, 2] = 0

					#--- this part is for gain of kmers
					
					for m in range(4): # iterate through all mutated nucleotides at this position
						old_idx = cn_convert(seq[pos])
						if m != old_idx: # not the original nucleotide 
							gain_kmers = set() # this should be evaluated for each mutated outcome (gain of kmers)

							for i in range(pos-args.l+1, pos+1): # iterate through all relevant positions
								if (i >= 0) and (i + args.l <= len(seq)): 
									original_kmer = seq[i:(i+args.l)]

									if re.search(r'[^ACGT]', original_kmer) is None: # make sure it only contains ACGT	
										nt_pos_in_kmer = pos - i
										new_kmer = original_kmer[:nt_pos_in_kmer] + cn_convert(m) + original_kmer[nt_pos_in_kmer+1:]
										gain_kmers.add(new_kmer)

							if len(gain_kmers) != 0 and (not gain_kmers & kmer_excl_gain): # none of kmers gained due to the specific mutation at this position is in the exclusion list

								for i in range(pos-args.l+1, pos+1): # iterate through all relevant positions
									if (i >= 0) and (i + args.l <= len(seq)): 
										original_kmer = seq[i:(i+args.l)]

										if re.search(r'[^ACGT]', original_kmer) is None: # make sure it only contains ACGT	
											nt_pos_in_kmer = pos - i
											new_kmer = original_kmer[:nt_pos_in_kmer] + cn_convert(m) + original_kmer[nt_pos_in_kmer+1:]

											# update gain of new_kmer with running average and distance 
											if new_kmer in gain_dtemp:
												gain_dtemp[new_kmer][old_idx, nt_pos_in_kmer, 0] += 1
												avg_last = gain_dtemp[new_kmer][old_idx, nt_pos_in_kmer, 1]
												gain_dtemp[new_kmer][old_idx, nt_pos_in_kmer, 1] = avg_last + (tl_diff[pos,m] - avg_last) / gain_dtemp[new_kmer][old_idx, nt_pos_in_kmer, 0]
												gain_dtemp[new_kmer][old_idx, nt_pos_in_kmer, 2] = gain_dtemp[new_kmer][old_idx, nt_pos_in_kmer, 2] + (tl_diff[pos,m] - gain_dtemp[new_kmer][old_idx, nt_pos_in_kmer, 1]) * (tl_diff[pos,m] - avg_last)
											else:
												gain_dtemp[new_kmer] = np.zeros((4, len(new_kmer), 3))
												gain_dtemp[new_kmer][~reverse_one_hot(new_kmer),:] = np.nan
												gain_dtemp[new_kmer][old_idx, nt_pos_in_kmer, 0] = 1
												gain_dtemp[new_kmer][old_idx, nt_pos_in_kmer, 1] = tl_diff[pos,m]
												gain_dtemp[new_kmer][old_idx, nt_pos_in_kmer, 2] = 0
						
	return(gain_dtemp, loss_dtemp)
		
def add_to_first_dict(d1, d2, multi_dim = False):
	# d1 must have all keys that d2 have
	# https://math.stackexchange.com/questions/2971315/how-do-i-combine-standard-deviations-of-two-groups
	for key in d2:
		if key in d1:
			if not multi_dim:
				count_last = d1[key][0]
				d1[key][0] = d1[key][0] + d2[key][0]
				avg_last = d1[key][1]
				d1[key][1] = (count_last * d1[key][1] + d2[key][0] * d2[key][1]) / (count_last + d2[key][0])
				d1[key][2] = d1[key][2] + d2[key][2] + count_last * d2[key][0] * (avg_last - d2[key][1]) * (avg_last - d2[key][1]) / (count_last + d2[key][0])
			else:
				mask = reverse_one_hot(key)
				count_last = np.array(d1[key][mask,0])
				d1[key][mask,0] = d1[key][mask,0] + d2[key][mask,0]
				avg_last = np.array(d1[key][mask,1])
				d1[key][mask,1] = (count_last * d1[key][mask,1] + d2[key][mask,0] * d2[key][mask,1]) / (count_last + d2[key][mask,0])
				d1[key][mask,2] = d1[key][mask,2] + d2[key][mask,2] + count_last * d2[key][mask,0] * (avg_last - d2[key][mask,1]) * (avg_last - d2[key][mask,1]) / (count_last + d2[key][mask,0])
		else: 
			sys.exit('No key ' + key + ' in dictionary when excuting function "add_to_first_dict"!\n')

def remove_pseudo_count(d, multi_dim = False):
	for k in d:
		if multi_dim:
			nt_mask = reverse_one_hot(k)
			zero_count_mask = d[k][:,:,0] > 1
			mask = nt_mask & zero_count_mask

			count_last = np.array(d[k][mask,0])
			d[k][mask,0] = count_last - 1
			avg_last = np.array(d[k][mask,1])
			d[k][mask,1] = (d[k][mask,1] * count_last - 0) / d[k][mask,0]
			d[k][mask,2] = d[k][mask,2] - 0 - count_last / d[k][mask,0] * avg_last * avg_last

			zero_count_not_self_nt_mask = nt_mask & (~zero_count_mask)
			d[k][zero_count_not_self_nt_mask,:] = np.zeros(3)
		else:
			if d[k][0] > 1:
				count_last = d[k][0]
				d[k][0] = d[k][0] - 1
				avg_last = d[k][1]
				d[k][1] = (d[k][1] * count_last - 0) / d[k][0]
				d[k][2] = d[k][2] - 0 - count_last / d[k][0] * avg_last * avg_last		
			else:
				d[k] = [0, 0, 0]

def reverse_one_hot(kmer):
    # generate a boolen 2D array, with reverse 1-hot encoding of the kmer 
    ary = np.full((4, len(kmer)), True, dtype=bool)
    for i in range(len(kmer)):
        ary[cn_convert(kmer[i]),i] = False
    return(ary)

def get_summary_from_kmer_array(kmer, ary):
	# 'ary' 3D numpy array (4 [ACGT channels] x l [len of kmer] x 3 [count of instances, avg TL, dis TL])
	# for kmer_loss arrays, the nucleotide channel mean the nucleotide it mutated to, so 3 out 4 channels have values in each position
	# for kmer_gain arrays, the nucleotide channel mean the nucleotide it mutated from, so 3 out 4 channels have values in each position
	mask = reverse_one_hot(kmer)
	n = np.sum(ary[mask,0])
	if n > 0:
		avg = np.sum(ary[mask,0] * ary[mask,1]) / n
	else:
		avg = 0
	if n > 1:
		ss = (np.sum(ary[mask,2]) + np.sum(ary[mask,0] * (ary[mask,1] - avg) ** 2)) / (n - 1)
	else:
		ss = np.nan

	return(n, avg, ss)

def one_sample_t_test(avg, n, std, ref = 0, one_tailed = True):
	if n < 2:
		return(1)
	else:
		t_statistic = (avg - ref) / (std / (n ** 0.5))
		df = n - 1
		if one_tailed:
			p_value = 1 - stats.t.cdf(abs(t_statistic), df)
		else:
			p_value = 2 * (1 - stats.t.cdf(abs(t_statistic), df))
		return(p_value)

#
#####################################################################################################################
###------the script runs from there-------------------------------------------------------
t_start = time() # timer start

# parse the input
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest = 'i', type = str, help = 'input file in h5 format with CNN ISM result')
parser.add_argument('-r', '--reference', dest = 'r', type = str, help = 'reference poly(A) site file')
parser.add_argument('-l', '--length', dest = 'l', type = int, default = 6, help = 'length of kmers to examine')
parser.add_argument('-o', '--output', dest = 'o', type = str, help = 'output folder')
parser.add_argument('-u', '--utr_length', dest = 'u', type = int, help = 'longest utr region to examine (from 3\' ends)')
parser.add_argument('-n', '--name', dest = 'n', type = str, help = 'output file name prefix')
args = parser.parse_args()  

# output folder
os.makedirs(args.o, exist_ok=True)

if args.n is None:
	idf = args.i.split('/')[-1].split('.')[0]
else:
	idf = args.n

if args.u:
	idf = idf + f"_utr_len_limit_{args.u}"

f_log = open(os.path.join(args.o, idf + '_kmer_iterative_exclusion_pA_weight_avg_dis_log.txt'), 'w')

# folder for storing output data of each iterative round
temp_dir = os.path.join(args.o, 'Iterative_exclusion_output_by_rounds')
os.makedirs(temp_dir, exist_ok=True)

pwrite(f_log, 'Input file: ' + args.i)
pwrite(f_log, 'Output data in this folder:\n' + args.o)
pwrite(f_log, '\nExamining kmers with length of ' + str(args.l) + '...' + timer())

# make a dictionary of poly(A) sites for filtering those are uncertain (low confidence)
if args.r:
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
	pwrite(f_log, 'Number of confident pA sites in the reference file: ' + str(len(pa_dict)) + timer())
else:
	pwrite(f_log, 'No reference file provided. Entries are not filtered.' + timer())

with h5py.File(args.i, 'r') as f, \
	open(os.path.join(args.o, idf + '_gain_kmer_iterative_exclusion_pA_weight_avg_dis.txt'), 'w') as f_out_gain, \
	open(os.path.join(args.o, idf + '_loss_kmer_iterative_exclusion_pA_weight_avg_dis.txt'), 'w') as f_out_loss:
	
	# get the poly(A) site ID
	if args.r:
		if type(f['id'][:][0]) == bytes:
			idx_sele_lst = [i for i in range(f['id'][:].shape[0]) if f['id'][:][i].decode('utf-8') in pa_dict]
		else:
			idx_sele_lst = [i for i in range(f['id'][:].shape[0]) if f['id'][:][i] in pa_dict]
	else:
		idx_sele_lst = f['id'][:]

	pwrite(f_log, '\nNumber of confident pA sites in the input file: ' + str(len(idx_sele_lst)) + timer())
	
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

	f_out_gain.write('\t'.join(['kmer', 'n', 'avg', 'dis', 'padj', 'round']) + '\n')
	f_out_loss.write('\t'.join(['kmer', 'n', 'avg', 'dis', 'padj', 'round']) + '\n')

	# make a set for kmers to exlcude
	kmer_excl_gain = set()
	kmer_excl_loss = set()
	kmer_gain_lst = []
	kmer_loss_lst = []
	kmer_gain_ary_lst = []
	kmer_loss_ary_lst = []
	n_round = 0
	flag_continue_gain = True
	flag_continue_loss = True
	while(True):
		if flag_continue_gain or flag_continue_loss:
			n_round += 1
			pwrite(f_log, '\n' + ('').join(['-']*100))
			pwrite(f_log, 'Round ' + str(n_round) + ' has started...' + timer())

			# make two master dictionary with kmer as keys, one is for gaining of a kmer and the other is for losing a kmer
			# {kmer: 3D numpy array (4 [ACGT channels] x l [len of kmer] x 3 [count of instances, avg TL, dis TL])}
			# "avg" is the mean and "dis" is std^2 * (n-1)
			# it's necessary to initiate with pseudo count at all position, otherwise dictionary combining in "add_to_first_dict" will have a "divide by zero error"
			dict_gain = {} 
			dict_loss = {}
			for kmer in generate_kmers(kmer_len_min = args.l, kmer_len_max = args.l):
				if kmer not in kmer_excl_gain:
					dict_gain[kmer] = np.zeros((4, len(kmer), 3))
					dict_gain[kmer][reverse_one_hot(kmer), 0] = 1
					dict_gain[kmer][~reverse_one_hot(kmer), :] = np.nan
				if kmer not in kmer_excl_loss:
					dict_loss[kmer] = np.zeros((4, len(kmer), 3))
					dict_loss[kmer][reverse_one_hot(kmer), 0] = 1
					dict_loss[kmer][~reverse_one_hot(kmer), :] = np.nan

			with concurrent.futures.ProcessPoolExecutor(n_threads) as pool:
				futures = pool.map(exam_kmers_mp, np.array_split(range(len(idx_sele_lst)), n_threads))
				for d_gain, d_loss in futures: # combine data from outputs from all processes
					add_to_first_dict(dict_gain, d_gain, multi_dim = True)
					add_to_first_dict(dict_loss, d_loss, multi_dim = True)			

			# need to remove the original pseudo count, which all have average weight 0 and standard deviation 0
			remove_pseudo_count(dict_gain, multi_dim = True)
			remove_pseudo_count(dict_loss, multi_dim = True)
			pwrite(f_log, 'Finished examining kmers...' + timer())	

			if flag_continue_gain:
				with open(os.path.join(temp_dir, idf + '_kmer_gain_round_' + str(n_round) + '.txt'), 'w') as f_gain_temp:
					summaries = {kmer: get_summary_from_kmer_array(kmer, ary) for kmer, ary in dict_gain.items()}
					sorted_kmers = sorted(dict_gain.items(), key=lambda x: summaries[x[0]][1], reverse=True)
					for i, (kmer, ary) in enumerate(sorted_kmers):
						lst = summaries[kmer]
						f_gain_temp.write(kmer + '\t' + '\t'.join(list(map(str, lst))) + '\n')
						if i == 0:
							if int(lst[0]) > 1:
								p_value = one_sample_t_test(avg = float(lst[1]), n = int(lst[0]), std = (float(lst[2]) / (int(lst[0])-1)) ** 0.5, ref = 0)
								padj = p_value * (len(dict_gain) - len(kmer_excl_gain))
								if padj > 0.05 or float(lst[1]) <= 0:
									flag_continue_gain = False
							else:
								flag_continue_gain = False
							
							if not flag_continue_gain:
								pwrite(f_log, '\nThe top kmer in the GAIN-of-kmer list did not pass the t-test for this round.')
							else:
								kmer_excl_gain.add(kmer)
								kmer_gain_lst.append(kmer)
								kmer_gain_ary_lst.append(ary)
								f_out_gain.write(kmer + '\t' + '\t'.join(list(map(str, lst))) + '\t' + str(padj) + '\t' + str(n_round) + '\n')
								pwrite(f_log, '\nThe top kmer in the GAIN-of-kmer list in this round is: ' + kmer + 
									' (avg = ' + f'{lst[1]:.3g}' + ', padj = ' + f'{padj:.3g}' + ')')
								pwrite(f_log,  pd.DataFrame(ary[:,:,1], index = list('ACGT'), columns = list(kmer)).to_string(float_format="%.3f"))
			else:
				pwrite(f_log, '\nNo kmers in the GAIN-of-kmer list have passed the t-test for the previous rounds.')


			if flag_continue_loss:				
				with open(os.path.join(temp_dir, idf + '_kmer_loss_round_' + str(n_round) + '.txt'), 'w') as f_loss_temp:
					summaries = {kmer: get_summary_from_kmer_array(kmer, ary) for kmer, ary in dict_loss.items()}
					sorted_kmers = sorted(dict_loss.items(), key=lambda x: summaries[x[0]][1])
					for i, (kmer, ary) in enumerate(sorted_kmers):
						lst = summaries[kmer]
						f_loss_temp.write(kmer + '\t' + '\t'.join(list(map(str, lst))) + '\n')
						if i == 0:
							if int(lst[0]) > 1:
								p_value = one_sample_t_test(avg = float(lst[1]), n = int(lst[0]), std = (float(lst[2]) / (int(lst[0])-1)) ** 0.5, ref = 0)
								padj = p_value * (len(dict_loss) - len(kmer_excl_loss))
								if padj > 0.05 or float(lst[1]) >= 0:
									flag_continue_loss = False
							else:
								flag_continue_loss = False

							if not flag_continue_loss:
								pwrite(f_log, '\nThe top kmer in the LOSS-of-kmer list did not pass the t-test for this round.')
							else:
								kmer_excl_loss.add(kmer)
								kmer_loss_lst.append(kmer)
								kmer_loss_ary_lst.append(ary)
								f_out_loss.write(kmer + '\t' + '\t'.join(list(map(str, lst))) + '\t' + str(padj) + '\t' + str(n_round) + '\n')
								pwrite(f_log, '\nThe top kmer in the LOSS-of-kmer list in this round is: ' + kmer + 
									' (avg = ' + f'{lst[1]:.3g}' + ', padj = ' + f'{padj:.3g}' + ')')
								pwrite(f_log,  pd.DataFrame(ary[:,:,1], index = list('ACGT'), columns = list(kmer)).to_string(float_format="%.3f"))
			else:
				pwrite(f_log, '\nNo kmers in the LOSS-of-kmer list have passed the t-test for the previous rounds.')
		else:
			pwrite(f_log, '\nNo kmers have passed the t-test in the previous round. The run will stop.')
			break

pwrite(f_log, 'Output data to h5 files...' + timer())
with h5py.File(os.path.join(args.o, idf + '_gain_kmer_iterative_exclusion_pA_weight_avg_dis.h5'), 'w') as f:
	f.create_dataset('kmer', (len(kmer_gain_lst),), dtype = h5py.special_dtype(vlen=str), compression = 'lzf')
	f['kmer'][:] = list(kmer_gain_lst)
	f.create_dataset('stats', (len(kmer_gain_ary_lst), 4, args.l, 3), compression = 'lzf')
	f['stats'][:] = np.stack(kmer_gain_ary_lst, axis = 0)

with h5py.File(os.path.join(args.o, idf + '_loss_kmer_iterative_exclusion_pA_weight_avg_dis.h5'), 'w') as f:
	f.create_dataset('kmer', (len(kmer_loss_lst),), dtype = h5py.special_dtype(vlen=str), compression = 'lzf')
	f['kmer'][:] = kmer_loss_lst
	f.create_dataset('stats', (len(kmer_loss_ary_lst), 4, args.l, 3), compression = 'lzf')
	f['stats'][:] = np.stack(kmer_loss_ary_lst, axis = 0)

pwrite(f_log, 'Final: ' + timer())		
f_log.close()
 
