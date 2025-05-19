import sys, concurrent.futures, math, argparse, re, os, matplotlib, h5py, yaml, optuna, gzip
import numpy as np
import pandas as pd
import tensorflow as tf

from time import time
from datetime import datetime
from tensorflow import keras
from sklearn.preprocessing import PowerTransformer
from scipy import stats
from collections import OrderedDict
from tensorflow.keras.models import load_model

matplotlib.use('Agg')

# Default global parameters can be updated by the input yaml file with the '-p' option
# some parameters can also be updated by the input options and they will override the default values and values provided in the yaml file
global_config = OrderedDict([
	('NAME', 'Global parameters'),
	('fa_string_filter', 'uncertain'), # if fasta entry contains this string, it will be filtered out
	('flag_input_header', True), # whether the input file with target values has a header
	('seq_3end_append', 'AAA'), # sequence appended at the 3'-end
	('len_max', 1000), # maximal length of sequence (from 3' ends, after appending the constant region, if applicable). This can also be updated with the '-l' option.
	('len_min', 10), # minimal length of 3' UTR sequence to be included, from 3' ends (those smaller than this won't be analyzed)
	('y_mode', 'none'), # how to transform target (tail length) value, can also be changed by input option '-t'
	('y_offset', 0), # an offset value applied to the target (tail length) value. Automatically modified based on the input data.
	('tl_cons', 35), # tail length constant to be subtracted if 'y_mode' is set to 'diff'
	('flag_minn', False), # boolean flag for whether the model is an MINN model 
	('predict_model_index', 1), # the index of group in the model for prediction (n-th model = index + 1)
	('batch_size', 128),
	('n_threads', 20),
	('chunk', 10000),
	('chunk_percent', 5)
])

def timer(t_start): 
	# a function to calculate runtime
	temp = str(time()-t_start).split('.')[0]
	temp =  '\t' + temp + 's passed...' + '\t' + str(datetime.now())
	return(temp)

def pwrite(f, text, timestamp=None):
	"""Write printed output to a file and print it to console."""
	log_text = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - {text}" if timestamp else text
	if f:
		f.write(log_text + '\n')
		f.flush() 
	print(log_text)

def make_log_file(filename, p_params = False, p_vars = False):
	# a function to generate log file
	f_log = open(filename, 'w')
	if isinstance(p_params, dict):
		pwrite(f_log, 'Default parameters:')
		for param in p_params:
			pwrite(f_log, f'\t{param}: {p_params[param]}')
	if isinstance(p_vars, dict):
		pwrite(f_log, '\nInput arguments:')
		for var in p_vars:
			pwrite(f_log, f'\t{var}: {p_vars[var]}')
	return(f_log)
	
def load_keras_models(path):
	if os.path.isdir(path):
		# If it's a folder, load all models in the folder (except the callback models)
		models = []
		for file_name in os.listdir(path):
			if 'callback' not in file_name:
				file_path = os.path.join(path, file_name)
				if os.path.isfile(file_path) and file_name.endswith('.keras'):
					models.append(load_model(file_path))
		return(models)

	elif os.path.isfile(path) and path.endswith('.keras'):
		# If it's a file, load the model
		return([load_model(path)])
	else:
		raise ValueError("The provided path is neither a valid folder nor a Keras model file.")

def fasta_id_parser(line, sep = '\t'):
	# a function to parse the id line in a fasta file
	lst = line.rstrip().lstrip('>').split(sep)
	return(lst)

def fasta_to_dict(file, string_filter = None, key_parser_sep = '\t', key_parser_pos = 0):
	# a function to make a dictionary from a fasta file 
	# use "string_filter" to filter out entries 
	# use "key_parser_sep" and "key_parser_pos" to determine how to extract the keys from the id lines
	temp_dict = {}
	temp_seq = ''
	current_id = None
	with open(file, 'r') as f:
		flag_seq_append = False
		while(True):
			line = f.readline()
			if line:
				if line[0] == '>':
					if current_id is not None:
						temp_dict[current_id] = temp_seq
					if (string_filter is None) or (string_filter not in line):
						current_id = line.rstrip().lstrip('>').split(key_parser_sep)[key_parser_pos]
						temp_seq = ''
						if current_id in temp_dict:
							sys.exit('Error! Input file (' + file + ') has non-unique IDs!')
					else:
						current_id = None
				else:
					if current_id is not None:
						temp_seq = temp_seq + line.strip()
			else:
				# last entry
				if current_id is not None:
					temp_dict[current_id] = temp_seq 
				break
	return(temp_dict)

def encode_seq(seq):
	# convert a DNA sequence to one-hot encoding matrix
	mat = np.zeros((len(seq), 4))
	for i in range(len(seq)):
		if seq[i] == 'A':
			mat[i,0] = 1.0
		elif seq[i] == 'C':
			mat[i,1] = 1.0
		elif seq[i] == 'G':
			mat[i,2] = 1.0
		elif seq[i] == 'T' or seq[i] == 'U':
			mat[i,3] = 1.0
	return(mat) 

def decode_seq(ary):
	# ary must be a 2D array, with first dimension along the seq and second dimension encoding nucleotide
	nt_lst = ['N', 'A', 'C', 'G', 'T']
	seq = ''
	for i in range(ary.shape[0]):
		seq = seq + nt_lst[int(np.dot(ary[i,0:4], np.arange(1,5)))]
	return(seq)

def mutate_sequence(seq_ary):
	'''
	mutate the matrix for every position and every possible nucleotide
	shape for 'seq_ary': length of the sequence * 5 (A, C, G, T, N channels)
	'''
	mut_seq_ary = []
	mut_pos_lsts = []
	for pos in range(seq_ary.shape[0]):
		if (seq_ary.shape[1] == 5 and seq_ary[pos,4] == 1) or seq_ary.shape[1] == 4: # must be 3' UTR
			if any(seq_ary[pos,:4] > 0): # cannot be 'N'
				original_idx = np.argmax(seq_ary[pos, :4])  # Get the original nucleotide
				mutated_idx = [idx for idx in range(4) if idx != original_idx]
				
				for idx in mutated_idx:
					mut_ary = np.copy(seq_ary)
					mut_ary[pos,:4] = 0
					mut_ary[pos, idx] = 1
					mut_seq_ary.append(mut_ary)
					mut_pos_lsts.append([pos, idx])

	return(np.stack(mut_seq_ary, axis = 0), mut_pos_lsts)

def process_seqs(line_lst):
	# convert sequences to features as X values and tail length to Y values
	X_mat_lst = []
	id_lst = []
	n_drop_fa_string_filter = 0 # number of entries dropped because the name of each sequence entry in the input file contains a string ("uncertain" annotations)
	n_drop_short_utr = 0 # number of entries dropped because of short UTR sequence
	n_drop_unknown_nt = 0 # numbner of entries dropped becasue of unknown characters in the sequence
	n_utr_w_cds = 0 # number of UTRs with CDS sequences appended
	global args

	for lines in line_lst:
		sele_id = fasta_id_parser(lines[0], sep = '__')[0]
		if len(fasta_id_parser(lines[0], sep = '__')) > 1:
			transcript_id = fasta_id_parser(lines[0], sep = '__')[1] 
		else:
			transcript_id = fasta_id_parser(lines[0], sep = '__')[0]

		if (global_config['fa_string_filter'] is None) or (global_config['fa_string_filter'] not in lines[0]): # exclude "uncertain" annotations
			seq = lines[1].rstrip()
			
			# append the 3' end sequence (default: 3 nt from the poly(A) tail)
			seq = seq + global_config['seq_3end_append']
			
			if len(seq) >= global_config['len_min']: # UTR must be longer than this

				# truncate long UTR sequences and also pad short UTR sequences with N at the 5' ends
				utr_seq = seq 
				if len(seq) >= global_config['len_max']: 
					seq = seq[-global_config['len_max']:] 
				else:	
					seq = 'N' * (global_config['len_max'] - len(seq)) + seq

				# add the CDS track 
				if args.c:
					if len(utr_seq) >= global_config['len_max']:
						anno_track = [1 for i in range(len(seq))]
					else:
						anno_track = [0 for i in range(global_config['len_max']-len(utr_seq))] + [1 for i in range(len(utr_seq))]
						
					# add the CDS sequence if provided
					if 'cds_dict' in globals():
						if transcript_id in cds_dict or len(cds_dict) == 1: # if only one entry in the CDS dictionary, treat this as universal CDS for all UTR3 variants
							n_utr_w_cds += 1
							if transcript_id in cds_dict:
								cds_seq = cds_dict[transcript_id]
							else:
								cds_seq = list(cds_dict.values())[0]
								
							if len(cds_seq) > global_config['len_max']-len(utr_seq):
								seq = (cds_seq + utr_seq)[-global_config['len_max']:]
							else:
								seq = 'N'*(global_config['len_max']-len(cds_seq) - len(utr_seq)) + cds_seq + utr_seq

				# encode the sequence
				seq = seq.upper().replace('U', 'T')
				if re.search('(^[ACGTNU]+$)',seq):
					seq_mat = encode_seq(seq)
				else:
					pwrite(f_log, 'Error! Unknown characters in the sequence of this entry: ' + sele_id)
					n_drop_unknown_nt += 1
					break

				# add the cds track if provided
				if args.c:
					seq_mat = np.hstack((seq_mat, np.asarray(anno_track).reshape((-1,1))))

				X_mat_lst.append(seq_mat)
				id_lst.append(sele_id)

			else:
				n_drop_short_utr += 1
		else:
			n_drop_fa_string_filter += 1
			
	if len(X_mat_lst) > 0:
		return(
			np.stack(X_mat_lst, axis = 0), 
			pd.DataFrame({'id':id_lst}), 
			n_drop_fa_string_filter, 
			n_drop_short_utr, 
			n_drop_unknown_nt,
			n_utr_w_cds)
	else:
		return(None)

######------------------------------------------------------------------------------------------------------------
t_start = time() # timer start
parser = argparse.ArgumentParser()
parser.add_argument('-u', '--utr', dest = 'u', type = str, help = 'Fasta file for 3\'-UTR sequences')
parser.add_argument('-c', '--cds', dest = 'c', nargs = '?', const = True, help = 'Fasta file for CDS sequences')
parser.add_argument('-m', '--model', dest = 'm', type = str, help = 'either a folder of the model files or a model file (keras format) for prediction')
parser.add_argument('-p', '--prefix', dest = 'p', type = str, help = 'output file prefix')
parser.add_argument('-t', '--transform', dest = 't', type = str, default = 'none', help = 'how to transform Y, can be "none", "diff", log", "sqrt"')
parser.add_argument('-v', '--verbose', dest = 'v', type = int, default = 0, help = 'verbose for CNN')
parser.add_argument('-o', '--outfolder', dest = 'o', type = str, default = './', help = 'folder for output files')
parser.add_argument('-l', '--length', dest = 'l', type = str, help = 'maximal length of sequence')
args = parser.parse_args()

###---------------------------------------------------------------------------------------------------------
# This part update parameters 

# make output folders
os.makedirs(args.o, exist_ok=True)
	
if args.u is None:
	parser.error('Input file with sequences is required!')
else:
	if args.p:
		idf = args.p
	else:
		idf = args.u.split('/')[-1].split('.')[0]

f_log = make_log_file(os.path.join(args.o, idf + '_complete_log.txt'), p_vars = vars(args))

# output updated parameters in the log file
pwrite(f_log, '\n' + ('').join(['-']*100))
if args.l:
	global_config['len_max'] = int(args.l)

pwrite(f_log, 'Global parameters:')
for p in global_config:
	pwrite(f_log, f'\t{p}: {global_config[p]}')

###---------------------------------------------------------------------------------------------------------
# load model and update 'len_max' in global_config before data preparation
# Set available GPUs
gpus = tf.config.list_physical_devices('GPU')
if gpus:
	pwrite(f_log, f'\nDetected {len(gpus)} GPU cores.')
	for gpu in gpus:
		tf.config.experimental.set_memory_growth(gpu, True)
else:
	pwrite(f_log, f'\nNo GPU core is detected. Use CPU for this session. ')
gpus = tf.config.list_physical_devices('GPU')

# Load the pre-trained model
pwrite(f_log, '\n' + ('').join(['-']*100))
pwrite(f_log, 'Input model folder/file:\n' + args.m)
models = load_keras_models(args.m) 
pwrite(f_log, f'Number of models loaded: {len(models)}\n')
models[0].summary(print_fn = lambda x: pwrite(f_log, x))

model_input_shape = models[0].inputs[0].shape
pwrite(f_log, f'Model\'s input shape: {model_input_shape}')
if model_input_shape[1] != global_config['len_max']:
	sys.exit('Length of the input sequence (' + str(global_config['len_max']) + ') is different from the model input (' + str(model_input_shape[1]) + ')!')

if len(models[0].inputs) > 1:
	global_config['flag_minn'] = True
	pwrite(f_log, f"The input model is an MINN model. Group {global_config['predict_model_index'] + 1} will be used for predictions")

###---------------------------------------------------------------------------------------------------------
# This part prepares data 

pwrite(f_log, '\n' + ('').join(['-']*100))
pwrite(f_log, 'Prepare data for analysis...')

# if '-c' is present but the the CDS sequence file is not provided, a CDS annotation track is added to the input sequences arrays.
# if the CDS sequence file is not provided, make a dictionary 
if args.c:
	if args.c is True:
		pwrite(f_log, '\nNo CDS sequences are provided, but the CDS annotation track will be added.' + timer(t_start))
	else:
		pwrite(f_log, '\nInput CDS file:\n' + args.c)
		pwrite(f_log, '\nMake a dictionary of CDS sequences...' + timer(t_start))
		cds_dict = fasta_to_dict(args.c, key_parser_sep = '__', key_parser_pos = 0)
		pwrite(f_log, 'Number of entries in the fasta file: ' + str(len(cds_dict)))
else:
	pwrite(f_log, '\nNo CDS sequences are provided. Use only 3\'-UTR sequences...' + timer(t_start))

# convert sequence 
pwrite(f_log, '\nRead in the fasta sequence file...')
with open(args.u, 'r') as f:
	line_lst = []
	counting = 0
	n_drop_fa_string_filter = 0
	n_drop_short_utr = 0
	n_drop_unknown_nt = 0
	n_utr_w_cds = 0
	X_mat_lst = []
	Y_df_lst = []
	temp_seq = ''
	while(True):
		line = f.readline()
		if not line:
			line_lst.append([current_id_line, temp_seq]) # the last entry
			with concurrent.futures.ProcessPoolExecutor(min(global_config['n_threads'], len(line_lst))) as pool:
				futures = pool.map(process_seqs, np.array_split(line_lst, min(global_config['n_threads'], len(line_lst))))
				for future in futures:
					if future is not None:
						ary, df, n1, n2, n3, n4 = future
						X_mat_lst.append(ary)
						Y_df_lst.append(df)
						n_drop_fa_string_filter += n1
						n_drop_short_utr += n2
						n_drop_unknown_nt += n3
						n_utr_w_cds += n4
				pwrite(f_log, str(counting) + ' sequences processed...' + timer(t_start))
			break
		else:
			if line[0] == '>':
				counting += 1
				if temp_seq != '':
					line_lst.append([current_id_line, temp_seq])
					temp_seq = ''
				current_id_line = line
				if counting % global_config['chunk'] == 0:
					with concurrent.futures.ProcessPoolExecutor(global_config['n_threads']) as pool:
						futures = pool.map(process_seqs, np.array_split(line_lst, global_config['n_threads']))
						for future in futures:
							if future is not None:
								ary, df, n1, n2, n3, n4 = future
								X_mat_lst.append(ary)
								Y_df_lst.append(df)
								n_drop_fa_string_filter += n1
								n_drop_short_utr += n2
								n_drop_unknown_nt += n3
								n_utr_w_cds += n4
						pwrite(f_log, str(counting) + ' sequences processed...' + timer(t_start))
					line_lst = []
			else:
				temp_seq += line.strip()

pwrite(f_log, 'Total number of sequences in the input file: ' + str(counting))
pwrite(f_log, 'Number of sequences removed because sequences in the fasta file is uncertain: ' + str(n_drop_fa_string_filter))
pwrite(f_log, 'Number of sequences removed because 3\'-UTR seqeunces shorter than ' + str(global_config['len_min']) +': ' + str(n_drop_short_utr))
pwrite(f_log, 'Number of sequences removed because of un-recogenized characters in the sequence: ' + str(n_drop_unknown_nt))
if args.c:
	pwrite(f_log, 'Number of 3\'-UTRs with CDS sequences added to the 5\' end: ' + str(n_utr_w_cds))

X_ary = np.vstack(X_mat_lst)
Y_df = pd.concat(Y_df_lst, ignore_index = True)
pwrite(f_log, 'Number of sequences remained after filtering and encoding: ' + str(X_ary.shape[0]))

# check dimensionality
pwrite(f_log, '\nData matrix dimension, X: ' + str(X_ary.shape) + ' Y: ' + str(Y_df.shape))
if Y_df.shape[0] != X_ary.shape[0]:
	pwrite(f_log, 'Warining! Dimension of input matrices and the number of target values don\'t match!')

###---------------------------------------------------------------------------------------------------------
# This section performs in silico mutagenesis

pwrite(f_log, '\n' + ('').join(['-']*100))
pwrite(f_log, 'Performing in silico mutagenesis...' + timer(t_start))
pwrite(f_log, 'Total number of sequences for mutagenesis: ' + str(X_ary.shape[0]))
n_total = X_ary.shape[0]
feed_counter = 1
out_pred = os.path.join(args.o, idf + '_in_silico_mutagenesis_prediction_results.hdf5')
ism_score_lst = []

with h5py.File(out_pred, 'w') as f:
	f.create_dataset('id', (n_total,), dtype = h5py.special_dtype(vlen=str), compression = 'lzf')
	f['id'][:] = Y_df['id'].to_numpy()

	# output sequences (converted back from X ary)
	f.create_dataset('seq', (n_total,), dtype = h5py.special_dtype(vlen=str), compression = 'lzf')
	f['seq'][:] = np.asarray([decode_seq(X_ary[x,:,:]) for x in range(X_ary.shape[0])])
	
	# output predicted tail lengths
	f.create_dataset('wt_scores',(n_total,), compression = 'lzf')
	pred_lst = []
	for model in models:
		if global_config['flag_minn']:
			pred_lst.append(model.predict([X_ary, np.repeat(global_config['predict_model_index'], X_ary.shape[0])], verbose = args.v, batch_size = global_config['batch_size']).ravel())
		else:
			pred_lst.append(model.predict(X_ary, verbose = args.v, batch_size = global_config['batch_size']).ravel())
	f['wt_scores'][:] = np.mean(np.stack(pred_lst, axis = 0), axis = 0)
	pwrite(f_log, f'\nFinished predicting {X_ary.shape[0]} wild-type sequences.' + timer(t_start))

	# output CDS or 3'-UTR labels for each nucleotide in each sequence 
	if X_ary.shape[2] == 5:
		f.create_dataset('utr_bool', (n_total, global_config['len_max']), dtype = 'i1', compression = 'lzf')
		f['utr_bool'][:] = X_ary[:,:,4]

	# output predicted tail length for all in silico mutagenesis 
	f.create_dataset('all_scores', (n_total, global_config['len_max'], 4), compression = 'lzf')
	for i in range(X_ary.shape[0]):
		# mutate each nucleotide of this sequence
		X_mut_ary, mut_pos_lsts = mutate_sequence(X_ary[i,:,:])

		# prediction with each model and then calculate the average
		pred_lst = []
		for model in models:
			if global_config['flag_minn']:
				pred_lst.append(model.predict([X_mut_ary, np.repeat(global_config['predict_model_index'], X_mut_ary.shape[0])], verbose = args.v, batch_size = global_config['batch_size']).ravel())
			else:
				pred_lst.append(model.predict(X_mut_ary, verbose = args.v, batch_size = global_config['batch_size']).ravel())
		predictions = np.mean(np.stack(pred_lst, axis = 0), axis = 0)

		# fill a all-NaN array with predictions
		ism_scores = np.full((X_ary.shape[1], 4), np.nan)
		for j in range(len(mut_pos_lsts)):
			pos, mut_idx = mut_pos_lsts[j]
			ism_scores[pos, mut_idx] = predictions[j]

		ism_score_lst.append(ism_scores)
		if i >= n_total * global_config['chunk_percent'] / 100.0 * feed_counter:
			pwrite(f_log, f"{feed_counter * global_config['chunk_percent']} percent sequences (n = {i+1}) have been processed..." + timer(t_start))
			feed_counter += 1
	pwrite(f_log, '100% percent sequences (n=' + str(i+1) + ') have been processed...' + timer(t_start))
	f['all_scores'][:] = np.stack(ism_score_lst, axis = 0)

pwrite(f_log, 'Finished...' + timer(t_start))
f_log.close()




