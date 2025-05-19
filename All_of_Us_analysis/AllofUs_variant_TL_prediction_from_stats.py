'''
v20240519
'''
import os, pysam, time, re, argparse, gzip, tarfile, sys
import numpy as np
import pandas as pd
import tensorflow as tf
import keras.backend as K

from time import time
from time import sleep
from datetime import datetime
from Bio.Seq import reverse_complement
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import load_model
from collections import OrderedDict

params = OrderedDict([
    #('flag_input_header', True), # whether the input file has a header line
    ('seq_append', 'AAA'), # sequence to append at the 3'-end of each sequence for prediction
    ('len_max', 2000), # maximal length of sequence (from 3' ends, after appending constant region, if applicable)
    ('len_min', 10), # minimal length of 3' UTR sequence to be included, from 3' ends (those smaller than this won't be analyzed)
    ('flag_cds_track', False), # boolean flag for whether a fifth CDS track is needed for model prediction
    ('flag_multi_heads', False), # boolean flag indicating whether the model is a multi-head model
    ('model_index', 1), # the index in a multi-head model
    ('batch_size', 128)
])

# some utility functions
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
    f_log = open(filename, 'w')
    if isinstance(p_params, dict):
        pwrite(f_log, 'Updated global parameters:')
        for param in p_params:
            pwrite(f_log, param + ': ' + str(p_params[param]))
    if isinstance(p_vars, dict):
        pwrite(f_log, '\nInput arguments:')
        for var in p_vars:
            pwrite(f_log, var + ': ' + str(p_vars[var]))
        pwrite(f_log, '\n') 
    return(f_log)

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

def get_variant_seq(line_lst, pa_seq_dict):
    '''
    get the variant 3'-UTR sequence 
    
        line_lst: a list of lines in the input file containing variant information
        pa_seq_dict: {pA_id : [sequences of exons of the wild type (reference allele) in right order]}
    '''
    global bed_dict, f_log
    variant_seq_lst = []
    variant_info_lsts = []
    for line in line_lst:
        v_chr, v_pos, v_ref, v_alt, v_af, v_ac, v_an, v_id, v_type, v_filter, pa, transcript_id = line.strip().split('\t')

        # iterate through all exons to find the exon where the variant is located
        exon_idx = -1
        for i in range(len(bed_dict[pa])):
            # Note that the v_pos is from a VCF file, thus it's 1-based. The annotation is a bed file, so it's 0-based.
            # Note that in cases when the ref allele is a long deletion, it can start outside the exon, leads to "v_pos < int(bed_dict[pa][i][1])"
            if int(bed_dict[pa][i][1]) <= (int(v_pos) - 1 + len(v_ref)) and (int(v_pos) - 1) < int(bed_dict[pa][i][2]):
                exon_idx = i
                exon_info = '_'.join([bed_dict[pa][i][j] for j in [0,1,2,5]])
                
                if bed_dict[pa][i][5] == '+':
                    ref_exon_seq = pa_seq_dict[pa][i]
                    alt_exon_seq = ref_exon_seq[:max((int(v_pos) - int(bed_dict[pa][i][1]) - 1),0)] + v_alt + ref_exon_seq[(int(v_pos) - int(bed_dict[pa][i][1]) - 1 + len(v_ref)):]
                else:
                    ref_exon_seq = reverse_complement(pa_seq_dict[pa][i])
                    alt_exon_seq = ref_exon_seq[:max((int(v_pos) - int(bed_dict[pa][i][1]) - 1),0)] + v_alt + ref_exon_seq[(int(v_pos) - int(bed_dict[pa][i][1]) - 1 + len(v_ref)):]
                    alt_exon_seq = reverse_complement(alt_exon_seq)

                break

        if exon_idx == -1:
            pwrite(f_log, f"This entry is not found in any annotated exons: {line}")
            sys.exit()

        # modify the exon to turn the wild type to the variant
        wt_seq = ''.join(pa_seq_dict[pa])
        alt_seq_lst = pa_seq_dict[pa][:]
        alt_seq_lst[exon_idx] = alt_exon_seq
        alt_seq = ''.join(alt_seq_lst)
        variant_seq_lst.append(alt_seq)

        # calculate the distance between the variant position and the 3'-end
        # note that this metric is meaningful only when the variant type is SNV
        if bed_dict[pa][exon_idx][5] == '+':
            dis2end = max(int(bed_dict[pa][exon_idx][2]) - int(v_pos) - len(v_ref) + 1, 0) + np.sum([len(seq) for n in range(len(alt_seq_lst)) if n > exon_idx])
        else:
            dis2end = int(v_pos) - int(bed_dict[pa][exon_idx][1]) - 1 + len(v_ref) - 1 + np.sum([len(seq) for n in range(len(alt_seq_lst)) if n > exon_idx])

        # check if the varaint has eliminated or created new CPE or PAS motifs
        n_cpe = len(get_motif_positions(alt_seq, 'TTTTA')) - len(get_motif_positions(wt_seq, 'TTTTA'))
        n_pas = len(get_motif_positions(alt_seq, 'AWTAAA')) - len(get_motif_positions(wt_seq, 'AWTAAA'))          

        # add information to the metadata dataframe
        variant_info_lsts.append([v_chr, v_pos, v_ref, v_alt, v_af, v_ac, v_an, v_id, v_type, dis2end, n_cpe, n_pas, bed_dict[pa][i][3], exon_info, pa])
        variant_df = pd.DataFrame(variant_info_lsts, columns = ['chromosome', 'position', 'ref', 'alt', 'af', 'ac', 'an', 'id', 'type', 'dis2end', 'cpe', 'pas', 'ref_pa', 'exon_info', 'pa_site'])

    return(variant_seq_lst, variant_df)

def convert_and_predict(seq_lst, pa_site_lst):
    # convert sequences to 1-hot encoding and predict tail length with the model
    global models, params
    assert len(seq_lst) == len(pa_site_lst)
    
    if len(seq_lst) > 0:
        X_lst = []
        for i, seq in enumerate(seq_lst):
            seq = seq + params['seq_append']
            if params['flag_cds_track']:
                utr_seq = seq
                
                global cds_dict, bed_dict
                transcript_id = bed_dict[pa_site_lst[i]][0][3].split('__')[1]
                if transcript_id in cds_dict:
                    cds_seq = cds_dict[transcript_id]
                else:
                    cds_seq = ''
                
                seq = cds_seq + utr_seq
                cds_track = [0] * len(cds_seq) + [1] * len(utr_seq)

            if len(seq) > params['len_max']: 
                    mat = encode_seq(seq[-params['len_max']:])
            else:
                mat = encode_seq('N'*(params['len_max'] - len(seq)) + seq)

            if params['flag_cds_track']:
                if len(cds_track) >= params['len_max']:
                    cds_track = np.array(cds_track[-params['len_max']:]).reshape((-1,1))
                else:
                    cds_track = np.array([0] * (params['len_max'] - len(cds_track)) + cds_track).reshape((-1,1))
                mat = np.hstack((mat, cds_track))

            X_lst.append(mat)

        X_mat = np.stack(X_lst, axis = 0)
        
        # predict with all models and take the average
        y_lst = []
        for model in models:
            if params['flag_multi_heads']:
                predictions = model.predict([X_mat, np.repeat(params['model_index'], X_mat.shape[0])], verbose = 0, batch_size = params['batch_size']).ravel()
            else:
                predictions = model.predict(X_mat, verbose = 0, batch_size = params['batch_size']).ravel()

            y_lst.append(predictions)

        y_avg = np.mean(np.stack(y_lst, axis = 0), axis = 0)
    return(y_avg)
 

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

######------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest = 'i', type = str, help = 'input files')
parser.add_argument('-b', '--bed', dest = 'b', type = str, help = 'bed file with 3\'-UTR annotations', \
    default = '/lab/solexa_bartel/coffee/Sequencing/Reference/Transcriptome/HS/Oocyte/UTR3/gencode.v25.annotation_annotation_fixed_UTR3_polished_by_PA_site_PolyA_DB.bed')
parser.add_argument('-f', '--fasta', dest = 'f', type = str, help = 'genome fasta file',\
    default = '/lab/solexa_bartel/coffee/Sequencing/Reference/Genome/HS/GRCh38.primary_assembly.genome.fasta')
parser.add_argument('-c', '--cds', dest = 'c', type = str, help = 'fasta file for the coding region')
parser.add_argument('-m', '--model', dest = 'm', type = str, help = 'either a folder of the model files or a model file (keras format) for prediction')
parser.add_argument('-n', '--name', dest = 'n', type = str, help = 'output file name prefix')
parser.add_argument('-o', '--out_folder', dest = 'o', type = str, default = './', help = 'output folder')
args = parser.parse_args()
t_start = time()

# construct the prefix if it is not in the input
if args.n:
    idf = args.n
else: 
    idf = args.i.split('/')[-1].split('.')[0]
os.makedirs(args.o, exist_ok=True)
f_log = make_log_file(os.path.join(args.o, idf + '_TL_INN_predictions.log'), p_params = params, p_vars = vars(args))

######------------------------------------------------------------------------------------------------------------
# Load the pre-trained model
pwrite(f_log, '\n' + ('').join(['-']*100))
pwrite(f_log, 'Predict with pre-trained models:')
pwrite(f_log, 'Input model folder/file:\n' + args.m)

try:
    models = load_keras_models(args.m) 
except:
    import inn_models
    models = load_keras_models(args.m) 

pwrite(f_log, f'Number of models loaded: {len(models)}\n')
models[0].summary(print_fn = lambda x: pwrite(f_log, x))

# check if the model is a simple model or a multi-head model
if len(models[0].inputs) > 1:
    params['flag_multi_heads'] = True
else:
    params['flag_multi_heads'] = False
pwrite(f_log, "params['flag_multi_heads'] has been updated: " + str(params['flag_multi_heads']))

# update 'len_max' based on model shape
model_input_shape = models[0].inputs[0].shape
pwrite(f_log, f'Model\'s input shape: {model_input_shape}\n')
if model_input_shape[1] != params['len_max']:
    pwrite(f_log, 'Length of the input sequence (' + str(params['len_max']) + ') is different from the model input (' + str(model_input_shape[1]) + ')!')
    params['len_max'] = model_input_shape[1]
    pwrite(f_log, "params['len_max'] has been updated: " + str(model_input_shape[1]))

# check if a CDS track is needed
if model_input_shape[2] == 5:
    params['flag_cds_track'] = True
    pwrite(f_log, "params['flag_cds_track'] has been updated: " + str(params['flag_cds_track']))

######------------------------------------------------------------------------------------------------------------
# get the 3'-utr annotation for each poly(A) site from a bed file
pwrite(f_log, '\n' + ('').join(['-']*100))

pwrite(f_log, 'Load the genome fasta file:\n\t' + args.f + timer(t_start))
genome_fa = pysam.FastaFile(args.f)

pwrite(f_log, '\nLoad the poly(A) sites bed file:\n\t' + args.b + timer(t_start))
bed_dict = {} # {poly(A) site id : list of annotation lines}
with open(args.b, 'r') as f:
    for line in f:
        lst = line.strip().split('\t')
        pa = lst[3].split('__')[0]
        if pa in bed_dict:
            bed_dict[pa].append(lst)
        else:
            bed_dict[pa] = [lst]
pwrite(f_log, 'Number of poly(A) sites: ' + str(len(bed_dict)))

######------------------------------------------------------------------------------------------------------------
# make a dictionary for the CDS if necessary
if params['flag_cds_track']:
    pwrite(f_log, '\n' + ('').join(['-']*100))
    pwrite(f_log, f"The model contains a track for the CDS annotation.")
    if args.c:
        pwrite(f_log, '\nInput CDS file:\n' + args.c)
        pwrite(f_log, '\nMake a dictionary of CDS sequences...' + timer(t_start))
        cds_dict = fasta_to_dict(args.c, key_parser_sep = '__', key_parser_pos = 0)
        pwrite(f_log, 'Number of entries in the fasta file: ' + str(len(cds_dict)))
    else:
        cds_dict = {}
        pwrite(f_log, '\nNo CDS sequences are provided. "N" nucleotides will be paded.' + timer(t_start))
    # cds_dict = {transcript_id : cds_seq}

######------------------------------------------------------------------------------------------------------------
# run through the variant file and obtain wild-type sequences
pwrite(f_log, '\n' + ('').join(['-']*100))
pwrite(f_log, 'Obtain wild-type sequences associated with variants...' + timer(t_start))
pa_seq_dict = {} # {pA_id : [sequences of exons of the wild type (reference allele) in right order]}
pa_site_lst = [] # [pA_id,...]
wt_seqs = [] # [sequence of the wild-type (reference allele),...]
with gzip.open(args.i, 'rt') as f:
    header = f.readline()
    for n_var, line in enumerate(f):
        pa = line.strip().split('\t')[-2] # the second to the last column is the poly(A) site in the input file
        if pa not in pa_seq_dict:
            # sort exons
            if (len(bed_dict[pa]) > 1):
                if bed_dict[pa][0][5] == '+':
                    bed_dict[pa] = sorted(bed_dict[pa], key = lambda x: int(x[1]))
                else:
                    bed_dict[pa] = sorted(bed_dict[pa], key = lambda x: int(x[1]), reverse = True)

            # get the reference 3'-UTR sequence
            wt_seq_lst = []
            for i in range(len(bed_dict[pa])):
                if bed_dict[pa][i][5] == '+':
                    seq = genome_fa.fetch(bed_dict[pa][i][0], int(bed_dict[pa][i][1]), int(bed_dict[pa][i][2])).upper()
                else:
                    seq = reverse_complement(genome_fa.fetch(bed_dict[pa][i][0], int(bed_dict[pa][i][1]), int(bed_dict[pa][i][2]))).upper()
                wt_seq_lst.append(seq)
            wt_seq = ''.join(wt_seq_lst)

            # the wild-type sequence must be at least this long
            if len(wt_seq + params['seq_append']) >= params['len_min']: 
                pa_seq_dict[pa] = wt_seq_lst
                pa_site_lst.append(pa)
                wt_seqs.append(wt_seq)

pwrite(f_log, 'Number of variants in the input file: ' + str(n_var + 1))
pwrite(f_log, 'Number of poly(A) sites in the input file: ' + str(len(pa_seq_dict)))
pwrite(f_log, 'Predict wild-type sequences...')
y_pred = convert_and_predict(wt_seqs, pa_site_lst)
wt_pred_df = pd.DataFrame({'pa_site':pa_site_lst, 'y_pred': y_pred})

######------------------------------------------------------------------------------------------------------------
# construct each variant sequence and use the model to predict the result
pwrite(f_log, '\n' + ('').join(['-']*100))
pwrite(f_log, 'Start predicting poly(A) tail length of the variants...' + timer(t_start))
out_df = pd.DataFrame(columns=['chromosome', 'position', 'ref', 'alt', 'af', 'ac', 'an', 'id', 'type', 'dis2end', 'cpe', 'pas', 'ref_pa', 'exon_info', 'pa_site', 'var_tl', 'ref_tl', 'diff_tl'])
out_lst = []
with gzip.open(args.i, 'rt') as f:
    header = f.readline()
    line_lst = []
    count = 0
    feed_counter = 1
    while(True):
        line = f.readline()
        if not line:
            if len(line_lst) > 0: # the last batch
                variant_seq_lst, variant_df = get_variant_seq(line_lst, pa_seq_dict)
                var_y_pred = convert_and_predict(variant_seq_lst, variant_df['pa_site'].to_numpy())
                variant_df['var_tl'] = var_y_pred
                variant_df['ref_tl'] = variant_df['pa_site'].map(wt_pred_df.set_index('pa_site')['y_pred'])
                variant_df['diff_tl'] = variant_df['var_tl'] - variant_df['ref_tl']
                out_lst.append(variant_df)
            pwrite(f_log, '100% percent sequences (n=' + str(count+1) + ') have been processed...' + timer(t_start))
            break
        else:
            lst = line.strip().split('\t')
            count += 1
            if lst[-2] in pa_seq_dict and lst[-3] == 'PASS': # variant must pass the filter
                line_lst.append(line)
                if len(line_lst) == params['batch_size']:
                    variant_seq_lst, variant_df = get_variant_seq(line_lst, pa_seq_dict)
                    var_y_pred = convert_and_predict(variant_seq_lst, variant_df['pa_site'].to_numpy())
                    variant_df['var_tl'] = var_y_pred
                    variant_df['ref_tl'] = variant_df['pa_site'].map(wt_pred_df.set_index('pa_site')['y_pred'])
                    variant_df['diff_tl'] = variant_df['var_tl'] - variant_df['ref_tl']
                    out_lst.append(variant_df)

                    # reset the line_lst
                    line_lst = []
            if count >= (n_var + 1) * 0.1 * feed_counter:
                pwrite(f_log, f"{feed_counter * 10} percent sequences (n = {count}) have been processed..." + timer(t_start))
                feed_counter += 1

out_df = pd.concat(out_lst, ignore_index = True)
out_df.to_csv(os.path.join(args.o, idf + '_INN_predictions.txt'), index = False, sep = '\t')
f_log.close()






