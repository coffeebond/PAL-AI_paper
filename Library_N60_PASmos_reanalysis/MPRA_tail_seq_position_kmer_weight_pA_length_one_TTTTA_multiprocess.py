# version 20240520
import sys, subprocess, numpy, time, math, concurrent.futures, argparse, itertools, re
import numpy as np
from time import time
from datetime import datetime

# parameters for the input file
header = True # if the input file has a header line 
col_seq = 1 # column number of the random sequence
col_weight = 2 # column number of the weight, in this case it's tail length

# parameters for counting kmers
seq_len = 60 # length of each sequence, for kmers at different positions

chunk = 1000000 # every number of lines to prompt the progress
chunk_lines = 50000 # number of lines to be process for each core
n_threads = 20 # number of cores to use

###------functions------------------------------------------------------------------------
def timer(): # calculate runtime
	temp = str(time()-t_start).split('.')[0]
	temp =  '\t' + temp + 's passed...' + '\t' + str(datetime.now())
	return temp

def generate_kmers(kmer_len):
    #Generates all kmers of length 'kmer_len' or shorter
    kmers = []
    for l in range(1, kmer_len+1):
        kmers += ["".join(kmer) for kmer in list(itertools.product(["A","C","G","T"], repeat=l))]
    return kmers

def pwrite(f, text):
	f.write(text + '\n')
	print(text)

def get_motif_positions(sequence, motif):
    # This function get positions of a motif in a query sequence (IUPAC code compatible)
    
    # Convert IUPAC code to regular expression
    iupac_regex = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "R": "[AG]",
        "Y": "[CT]",
        "S": "[GC]",
        "W": "[AT]",
        "K": "[GT]",
        "M": "[AC]",
        "B": "[CGT]",
        "D": "[AGT]",
        "H": "[ACT]",
        "V": "[ACG]",
        "N": "[ACGT]"
    }
    # Convert the motif to regular expression
    regex = "".join([iupac_regex.get(base, base) for base in motif.upper()])
    
    # Use regular expression to count motif occurrences in sequence
    matches = re.finditer(regex, sequence, re.IGNORECASE)
    
    # Use regular expression to count motif occurrences in sequence
    #count = len(re.findall(regex, sequence, re.IGNORECASE))
        
    return([match.start() for match in matches])

def exam_kmers_mp(lines):
	temp_dict = {}
	t2_dict = {}
	for line in lines:
		seq = line.rstrip().split('\t')[col_seq-1]

		# only examine variant that has one CPE
		if len(get_motif_positions(seq, 'TTTTA')) == 1:
			val = float(line.rstrip().split('\t')[col_weight-1])
			log_val = math.log(max(val, 1))
			
			for pos in range(len(seq)):
				for l in range(1, min((args.l+1), (len(seq)-pos + 1))):
					real_pos = pos + seq_len - len(seq) # real position (0-based) in the random region when aligned at the 3' end 
					real_pos_cpe = get_motif_positions(seq, 'TTTTA')[0] + seq_len - len(seq)
					rel_pos = pos - get_motif_positions(seq, 'TTTTA')[0]

					if rel_pos <= (-l) or rel_pos >= len('TTTTA'): # not overlapping with the CPE
						kmer = seq[pos:pos+l]
						if kmer in temp_dict:
							temp_dict[kmer][real_pos_cpe,real_pos,0] = temp_dict[kmer][real_pos_cpe,real_pos,0] + 1

							# update running average and distance
							avg_last = temp_dict[kmer][real_pos_cpe,real_pos,1]
							temp_dict[kmer][real_pos_cpe,real_pos,1] = avg_last + (val - avg_last) / temp_dict[kmer][real_pos_cpe,real_pos,0]
							temp_dict[kmer][real_pos_cpe,real_pos,2] = temp_dict[kmer][real_pos_cpe,real_pos,2] + (val - temp_dict[kmer][real_pos_cpe,real_pos,1]) * (val - avg_last)		
							log_avg_last = temp_dict[kmer][real_pos_cpe,real_pos,3]
							temp_dict[kmer][real_pos_cpe,real_pos,3] = log_avg_last + (log_val - log_avg_last) / temp_dict[kmer][real_pos_cpe,real_pos,0]
							temp_dict[kmer][real_pos_cpe,real_pos,4] = temp_dict[kmer][real_pos_cpe,real_pos,4] + (log_val - temp_dict[kmer][real_pos_cpe,real_pos,3]) * (log_val - log_avg_last)	
						else:
							temp_dict[kmer] = np.zeros((seq_len, seq_len, 5))
							temp_dict[kmer][real_pos_cpe,real_pos,0] = 1 
							temp_dict[kmer][real_pos_cpe,real_pos,1] = val 
							temp_dict[kmer][real_pos_cpe,real_pos,2] = 0 
							temp_dict[kmer][real_pos_cpe,real_pos,3] = log_val 
							temp_dict[kmer][real_pos_cpe,real_pos,4] = 0 

						idx = rel_pos+pos_idx_shift
						if kmer in t2_dict:
							t2_dict[kmer][idx,0] = t2_dict[kmer][idx,0] + 1

							# update running average and distance
							avg_last = t2_dict[kmer][idx,1]
							t2_dict[kmer][idx,1] = avg_last + (val - avg_last) / t2_dict[kmer][idx,0]
							t2_dict[kmer][idx,2] = t2_dict[kmer][idx,2] + (val - t2_dict[kmer][idx,1]) * (val - avg_last)		
							log_avg_last = t2_dict[kmer][idx,3]
							t2_dict[kmer][idx,3] = log_avg_last + (log_val - log_avg_last) / t2_dict[kmer][idx,0]
							t2_dict[kmer][idx,4] = t2_dict[kmer][idx,4] + (log_val - t2_dict[kmer][idx,3]) * (log_val - log_avg_last)	
						else:
							t2_dict[kmer] = np.zeros((2*seq_len+1, 5))
							t2_dict[kmer][idx,0] = 1 
							t2_dict[kmer][idx,1] = val 
							t2_dict[kmer][idx,2] = 0 
							t2_dict[kmer][idx,3] = log_val 
							t2_dict[kmer][idx,4] = 0 

	return(temp_dict, t2_dict)
		
def add_to_first_dict_2d(d1, d2):
	# d1 must have all keys that d2 have
	# https://math.stackexchange.com/questions/2971315/how-do-i-combine-standard-deviations-of-two-groups
	for key in d2:
		if key in d1:
			count_last = np.array(d1[key][:,0]) # be very careful with array slicing, it's different from list slicing
			d1[key][:,0] = d1[key][:,0] + d2[key][:,0]
			avg_last = np.array(d1[key][:,1])
			d1[key][:,1] = (count_last * d1[key][:,1] + d2[key][:,0] * d2[key][:,1]) / (count_last + d2[key][:,0])
			d1[key][:,2] = d1[key][:,2] + d2[key][:,2] + count_last * d2[key][:,0] * (avg_last - d2[key][:,1]) * (avg_last - d2[key][:,1]) / (count_last + d2[key][:,0])
			log_avg_last = np.array(d1[key][:,3])
			d1[key][:,3] = (count_last * d1[key][:,3] + d2[key][:,0] * d2[key][:,3]) / (count_last + d2[key][:,0])
			d1[key][:,4] = d1[key][:,4] + d2[key][:,4] + count_last * d2[key][:,0] * (log_avg_last - d2[key][:,3]) * (log_avg_last - d2[key][:,3]) / (count_last + d2[key][:,0])
		else: 
			sys.exit('No key ' + key + 'in dictionary when excuting function "add_to_first_dict_2d"!\n')

def add_to_first_dict_3d(d1, d2):
	# d1 must have all keys that d2 have
	# https://math.stackexchange.com/questions/2971315/how-do-i-combine-standard-deviations-of-two-groups
	for key in d2:
		if key in d1:
			count_last = np.array(d1[key][:,:,0]) # be very careful with array slicing, it's different from list slicing
			d1[key][:,:,0] = d1[key][:,:,0] + d2[key][:,:,0]
			avg_last = np.array(d1[key][:,:,1])
			d1[key][:,:,1] = (count_last * d1[key][:,:,1] + d2[key][:,:,0] * d2[key][:,:,1]) / (count_last + d2[key][:,:,0])
			d1[key][:,:,2] = d1[key][:,:,2] + d2[key][:,:,2] + count_last * d2[key][:,:,0] * (avg_last - d2[key][:,:,1]) * (avg_last - d2[key][:,:,1]) / (count_last + d2[key][:,:,0])
			log_avg_last = np.array(d1[key][:,:,3])
			d1[key][:,:,3] = (count_last * d1[key][:,:,3] + d2[key][:,:,0] * d2[key][:,:,3]) / (count_last + d2[key][:,:,0])
			d1[key][:,:,4] = d1[key][:,:,4] + d2[key][:,:,4] + count_last * d2[key][:,:,0] * (log_avg_last - d2[key][:,:,3]) * (log_avg_last - d2[key][:,:,3]) / (count_last + d2[key][:,:,0])
		else: 
			sys.exit('No key ' + key + 'in dictionary when excuting function "add_to_first_dict_3d"!\n')

#####################################################################################################################
###------the script runs from there-------------------------------------------------------
t_start = time() # timer start

# parse the input
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest = 'i', type = str, help = 'input file with pA lengths of all variants')
parser.add_argument('-l', '--length', dest = 'l', type = int, default = 5, help = 'longest kmer to examine')
args = parser.parse_args()  
f_out_1 = open(args.i.split('/')[-1].split('.txt')[0] + '_single_CPE_kmer_pos_pA_weight_avg_dis.txt', 'w')
f_out_2 = open(args.i.split('/')[-1].split('.txt')[0] + '_single_CPE_kmer_rel_pos_pA_weight_avg_dis.txt', 'w')
f_log = open(args.i.split('/')[-1].split('.txt')[0] + '_single_CPE_kmer_pos_pA_weight_avg_dis_log.txt', 'w')

pwrite(f_log, 'Examining position kmers up to length ' + str(args.l) + ' with respect to the sequence length of ' + str(seq_len) + '...' + timer())

# make a master dictionary with kmer as keys
# dimension 1: CPE positions
# dimension 2: kmer positions
# dimension 3: values for each position, including: count of sequences, avg TL, dis TL, avg log TL, dis log TL
# "avg" is the mean and "dis" is std^2 * (n-1)
# it's necessary to initiate with pseudo count at all position, otherwise dictionary combining in "add_to_first_dict" will have a "divide by zero error"
mdict = {} 
for kmer in generate_kmers(args.l):
	mdict.setdefault(kmer, np.zeros((seq_len, seq_len, 5)))
	mdict[kmer][:,:,0] = mdict[kmer][:,:,0] + 1

# make a second master dictionary with kmer as keys
# dimension 1: CPE-kmer relative positions
# dimension 2: values for each position, including: count of sequences, avg TL, dis TL, avg log TL, dis log TL
# note that there is a shift between relative position and the index in the 2nd dimension
# "avg" is the mean and "dis" is std^2 * (n-1)
# it's necessary to initiate with pseudo count at all position, otherwise dictionary combining in "add_to_first_dict" will have a "divide by zero error"
ndict = {} 
for kmer in generate_kmers(args.l):
	ndict.setdefault(kmer, np.zeros((seq_len * 2 + 1, 5)))
	ndict[kmer][:,0] = ndict[kmer][:,0] + 1
pos_idx_shift = 60 # index 0 is relative position -60


if args.i.endswith('.gz'):
	proc1 = subprocess.Popen(['zless', args.i], stdout = subprocess.PIPE)
elif args.i.endswith('.txt'):
	proc1 = subprocess.Popen(['less', args.i], stdout = subprocess.PIPE)
else:
	sys.exit('Wrong file type: ' + args.i)
proc2 = subprocess.Popen(['cut', '-f' + str(col_seq) + ',' + str(col_weight)], stdin = proc1.stdout, stdout = subprocess.PIPE)

if header:
	proc2.stdout.readline() # get rid of the header

counting = 0
rounds = 0
counting_sum = 0
chunk_temp = chunk
line_lst = []

while(True):
	line = proc2.stdout.readline().decode('utf-8')
	if not line:
		if len(line_lst) > 0 :
			with concurrent.futures.ProcessPoolExecutor(min(n_threads, len(line_lst))) as pool:
				futures = pool.map(exam_kmers_mp, np.array_split(line_lst, min(n_threads, len(line_lst))))
				for d3, d2 in futures: # combine data from outputs from all processes
					add_to_first_dict_3d(mdict, d3)
					add_to_first_dict_2d(ndict, d2)
				counting_sum += counting
				pwrite(f_log, str(counting_sum) + ' reads processed...' + timer())
		break
	else:
		counting += 1
		line_lst.append(line)
		if counting % (chunk_lines * n_threads) == 0:
			rounds += 1
			with concurrent.futures.ProcessPoolExecutor(n_threads) as pool:
				futures = pool.map(exam_kmers_mp, np.array_split(line_lst, n_threads))
				for d3, d2 in futures: # combine data from outputs from all processes
					add_to_first_dict_3d(mdict, d3)
					add_to_first_dict_2d(ndict, d2)
				counting_sum = counting * rounds
				if counting_sum >= chunk_temp:
					chunk_temp += chunk
					pwrite(f_log, str(counting_sum) + ' reads processed...' + timer())
			line_lst = []
			counting = 0

# need to remove the original pseudo counts, which all have average weight 0 and standard deviation 0
# if some positions still have only the pseudo count 1, reset them back to 0 in the end
f_out_1.write('\t'.join(['kmer', 'cpe_pos', 'kmer_pos', 'n', 'avg', 'dis', 'log_avg', 'log_dis']) + '\n')
for kmer in sorted(mdict.keys(), key = lambda x: (len(x), x)):
	count_last = np.array(mdict[kmer][:,:, 0])
	mdict[kmer][:,:, 0] = np.where(count_last > 1,  count_last - 1, count_last)
	avg_last = np.array(mdict[kmer][:,:, 1])
	mdict[kmer][:,:, 1] = (mdict[kmer][:,:, 1] * count_last - 0) / mdict[kmer][:,:, 0]
	mdict[kmer][:,:, 2] = mdict[kmer][:,:, 2] - 0 - count_last / mdict[kmer][:,:, 0] * avg_last * avg_last
	log_avg_last = np.array(mdict[kmer][:,:, 3])
	mdict[kmer][:,:, 3] = (mdict[kmer][:,:, 3] * count_last - 0) / mdict[kmer][:,:, 0]
	mdict[kmer][:,:, 4] = mdict[kmer][:,:, 4] - 0 - count_last / mdict[kmer][:,:, 0] * log_avg_last * log_avg_last
	mdict[kmer][:,:, 0] = np.where(count_last == 1,  0, mdict[kmer][:,:, 0])
	for cpe_pos in range(mdict[kmer].shape[0]):
		for kmer_pos in range(mdict[kmer].shape[1]):
			f_out_1.write(kmer + '\t' + str(cpe_pos+1) + '\t' + str(kmer_pos+1) + '\t' + '\t'.join(list(map(str, mdict[kmer][cpe_pos,kmer_pos,:]))) + '\n')

f_out_2.write('\t'.join(['kmer', 'rel_pos', 'n', 'avg', 'dis', 'log_avg', 'log_dis']) + '\n')
for kmer in sorted(ndict.keys(), key = lambda x: (len(x), x)):
	count_last = np.array(ndict[kmer][:, 0])
	ndict[kmer][:, 0] = np.where(count_last > 1,  count_last - 1, count_last)
	avg_last = np.array(ndict[kmer][:, 1])
	ndict[kmer][:, 1] = (ndict[kmer][:, 1] * count_last - 0) / ndict[kmer][:, 0]
	ndict[kmer][:, 2] = ndict[kmer][:, 2] - 0 - count_last / ndict[kmer][:, 0] * avg_last * avg_last
	log_avg_last = np.array(ndict[kmer][:, 3])
	ndict[kmer][:, 3] = (ndict[kmer][:, 3] * count_last - 0) / ndict[kmer][:, 0]
	ndict[kmer][:, 4] = ndict[kmer][:, 4] - 0 - count_last / ndict[kmer][:, 0] * log_avg_last * log_avg_last
	ndict[kmer][:, 0] = np.where(count_last == 1,  0, ndict[kmer][:, 0])
	for idx in range(ndict[kmer].shape[0]):
		rel_pos = idx - pos_idx_shift
		f_out_2.write(kmer + '\t' + str(rel_pos) + '\t' + '\t'.join(list(map(str, ndict[kmer][idx,:]))) + '\n')

pwrite(f_log, 'Final: ' + timer())		

 
