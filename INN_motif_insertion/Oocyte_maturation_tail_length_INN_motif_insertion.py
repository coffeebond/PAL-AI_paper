import sys, time, concurrent.futures, argparse, re, os, h5py
import numpy as np
import pandas as pd
import tensorflow as tf

from time import time
from datetime import datetime
from collections import OrderedDict
from tensorflow.keras.models import load_model

global_config = OrderedDict([
	('fa_string_filter', 'uncertain'), # if fasta entry contains this string, it will be filtered out
	('flag_input_header', True), # whether the input file with target values has a header
	('seq_3end_append', 'AAA'), # sequence appended at the 3'-end
	('len_max', 1056), #1056 # maximal length of sequence (from 3' ends, after appending constant region, if applicable)
	('len_min', 10), # minimal length of 3' UTR sequence to be included, from 3' ends (those smaller than this won't be analyzed)
	('tl_cons', 35), # tail length constant to be subtracted if 'y_mode' is set to 'diff'
	('flag_minn', False), # boolean flag for whether the model is an MINN model 
	('predict_model_index', 1), # the index of group in the model for prediction (n-th model = index + 1)
	('batch_size', 128),
	('n_threads', 4),
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
	
def motif_insert(seq_ary, motif):
	# insert the motif at every positions
	# shape of "seq_ary": (length of sequence, nucleotide identitiy with utr/cds track)

	motif_ary = encode_seq(motif)
	seq_length = seq_ary.shape[0]
	motif_len = len(motif)

	valid_positions = [
		pos for pos in range(seq_length - motif_len + 1)
		if (seq_ary.shape[1] == 4 or seq_ary[pos, 4] == 1)  # Ensure it's in 3' UTR
		and all(np.any(seq_ary[pos + i, :4] > 0) for i in range(motif_len))  # No 'N' bases
	]

	if not valid_positions:
		return np.array([]), []

	mut_seq_ary = np.tile(seq_ary, (len(valid_positions), 1, 1))

	for i, pos in enumerate(valid_positions):
		mut_seq_ary[i, pos:pos + motif_len, :4] = motif_ary  # Replace with motif

	return(mut_seq_ary, valid_positions)

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
parser.add_argument('-c', '--cds', dest = 'c', type = str, help = 'Fasta file for CDS sequences')
parser.add_argument('-f', '--motif', dest = 'f', type = str, help = 'Text file with motifs to insert')
parser.add_argument('-m', '--model', dest = 'm', type = str, help = 'CNN model file in h5 format')
parser.add_argument('-l', '--lenMax', dest = 'l', type = int, help = 'Maximal length of sequence (from 3\'-end) used in the analysis')
parser.add_argument('-o', '--output', dest = 'o', type = str, help = 'folder for output files')
parser.add_argument('-p', '--prefix', dest = 'p', type = str, help = 'output file prefix')
parser.add_argument('-v', '--verbose', dest = 'v', type = int, default = 0, help = 'verbose for CNN')
args = parser.parse_args()


###---------------------------------------------------------------------------------------------------------
# process input parameters
# output folder
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

##---------------------------------------------------------------------------------------------------------
# This part prepares data 

pwrite(f_log, '\n' + ('').join(['-']*100))
pwrite(f_log, 'Prepare data for analysis...')

###---------------------------------------------------------------------------------------------------------
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
# This section performs motif insertion

pwrite(f_log, '\n' + ('').join(['-']*100))
pwrite(f_log, '\nPerforming motif insertion analysis...' + timer(t_start))
pwrite(f_log, 'Total number of sequences for mutagenesis: ' + str(X_ary.shape[0]))

motif_lst = []
with open(args.f, 'r') as f:
	for m in f:
		motif_lst.append(m.strip())
pwrite(f_log, 'Motifs to be inserted: ')
pwrite(f_log, ', '.join(motif_lst))

n_total = X_ary.shape[0]

# wild-type predictions
pred_lst = []
for model in models:
	if global_config['flag_minn']:
		pred_lst.append(model.predict([X_ary, np.repeat(global_config['predict_model_index'], n_total)], verbose = args.v, batch_size = global_config['batch_size']).ravel())
	else:
		pred_lst.append(model.predict(X_ary, verbose = args.v, batch_size = global_config['batch_size']).ravel())
wt_predictions = np.mean(np.stack(pred_lst, axis = 0), axis = 0)
pwrite(f_log, f'\nFinished predicting {X_ary.shape[0]} wild-type sequences.' + timer(t_start))

for motif in motif_lst:
	pwrite(f_log, '\nInserting motif: ' + motif)
	out_pred = f"{args.o}{idf}_motif_insertion_{motif}_prediction_results.hdf5"

	with h5py.File(out_pred, 'w') as f:
		pwrite(f_log, 'Saving information...' + timer(t_start))
		# output id
		f.create_dataset('id', (n_total,), dtype = h5py.special_dtype(vlen=str), compression = 'lzf')
		f['id'][:] = Y_df['id'].to_numpy()
			
		# output sequences (converted back from X ary)
		f.create_dataset('seq', (n_total,), dtype = h5py.special_dtype(vlen=str), compression = 'lzf')
		f['seq'][:] = np.asarray([decode_seq(X_ary[x,:,:]) for x in range(n_total)])
		
		# output predicted tail lengths
		f.create_dataset('wt_scores',(n_total,), compression = 'lzf')
		f['wt_scores'][:] = wt_predictions
		
		# output CDS or 3'-UTR labels for each nucleotide in each sequence 
		if X_ary.shape[2] == 5:
			f.create_dataset('utr_bool', (n_total, global_config['len_max']), dtype = 'i1', compression = 'lzf')
			f['utr_bool'][:] = X_ary[:,:,4]

		# output predicted tail length for all in silico mutagenesis 
		grp = f.create_group('mi_scores')
		feed_counter = 1

		mi_score_lst = []
		X_mi_ary_lst = []
		motif_pos_lsts = []
		pwrite(f_log, 'Start motif insertions...' + timer(t_start))
		for i in range(X_ary.shape[0]):
			# for each sequence, get all matrices, predict Y and reshape into one 2D array
			mi_ary, motif_pos_lst = motif_insert(X_ary[i,:,:], motif)

			X_mi_ary_lst.append(mi_ary)
			motif_pos_lsts.append(motif_pos_lst)

			if np.sum([len(lst) for lst in motif_pos_lsts]) >= 1000 * global_config['batch_size']:
				X_mi_ary = np.vstack(X_mi_ary_lst)

				# prediction with each model and then calculate the average
				pred_lst = []
				for model in models:
					if global_config['flag_minn']:
						pred_lst.append(model.predict([X_mi_ary, np.repeat(global_config['predict_model_index'], X_mi_ary.shape[0])], verbose = args.v, batch_size = global_config['batch_size']).ravel())
					else:
						pred_lst.append(model.predict(X_mi_ary, verbose = args.v, batch_size = global_config['batch_size']).ravel())
				predictions = np.mean(np.stack(pred_lst, axis = 0), axis = 0)

				# fill a all-NaN array with predictions
				mi_scores = np.full((len(motif_pos_lsts), X_ary.shape[1] - len(motif) + 1), np.nan)
				pred_idx_start = 0
				for j in range(len(motif_pos_lsts)):
					mi_scores[j, np.array(motif_pos_lsts[j])] = predictions[pred_idx_start:(pred_idx_start + len(motif_pos_lsts[j]))]
					pred_idx_start = pred_idx_start + len(motif_pos_lsts[j])

				mi_score_lst.append(mi_scores)

				# reset lists
				X_mi_ary_lst = []
				motif_pos_lsts = []


			if i >= n_total * global_config['chunk_percent'] / 100.0 * feed_counter:
				pwrite(f_log, str(feed_counter * global_config['chunk_percent']) + ' percent data (n=' + str(i+1) + ') have been processed...' + timer(t_start))
				feed_counter += 1
		
		# the last batch
		X_mi_ary = np.vstack(X_mi_ary_lst)

		# prediction with each model and then calculate the average
		pred_lst = []
		for model in models:
			if global_config['flag_minn']:
				pred_lst.append(model.predict([X_mi_ary, np.repeat(global_config['predict_model_index'], X_mi_ary.shape[0])], verbose = args.v, batch_size = global_config['batch_size']).ravel())
			else:
				pred_lst.append(model.predict(X_mi_ary, verbose = args.v, batch_size = global_config['batch_size']).ravel())
		predictions = np.mean(np.stack(pred_lst, axis = 0), axis = 0)

		# fill a all-NaN array with predictions
		mi_scores = np.full((len(motif_pos_lsts), X_ary.shape[1] - len(motif) + 1), np.nan)
		pred_idx_start = 0
		for j in range(len(motif_pos_lsts)):
			mi_scores[j, np.array(motif_pos_lsts[j])] = predictions[pred_idx_start:(pred_idx_start + len(motif_pos_lsts[j]))]
			pred_idx_start = pred_idx_start + len(motif_pos_lsts[j])

		mi_score_lst.append(mi_scores)

		pwrite(f_log, '100% percent data (n=' + str(i+1) + ') have been processed...' + timer(t_start))
		
		grp.create_dataset(motif, (n_total, global_config['len_max'] - len(motif) + 1), compression = 'lzf')
		grp[motif][:] = np.vstack(mi_score_lst)
		feed_counter = 1

	pwrite(f_log, '\nFinished inserting motif: ' + motif + timer(t_start))

pwrite(f_log, 'Finished...' + timer(t_start))
f_log.close()

