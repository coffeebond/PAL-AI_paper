'''
v20240509
'''
import os, pysam, time, re, argparse
import numpy as np
import pandas as pd
import tensorflow as tf
import keras.backend as K

from time import time
from time import sleep
from datetime import datetime
from Bio.Seq import reverse_complement
from pysam import VariantFile
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import load_model
from collections import OrderedDict

params = OrderedDict([
    ('chr_sele', 'chrY'), # selected chromosome for analysis
    ('seq_append', 'AAA'), # sequence to append at the 3'-end of each sequence for prediction
    ('len_max', 2000), # 1056 # maximal length of sequence (from 3' ends, after appending constant region, if applicable)
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
parser.add_argument('-i', '--input', dest = 'i', type = str, help = 'input vcf files')
parser.add_argument('-b', '--bed', dest = 'b', type = str, help = 'bed file with 3\'-UTR annotations', \
    default = '/lab/solexa_bartel/coffee/Sequencing/Reference/Transcriptome/HS/Oocyte/UTR3/gencode.v25.annotation_annotation_fixed_UTR3_polished_by_PA_site_PolyA_DB.bed')
parser.add_argument('-f', '--fasta', dest = 'f', type = str, help = 'genome fasta file',\
    default = '/lab/solexa_bartel/coffee/Sequencing/Reference/Genome/HS/GRCh38.primary_assembly.genome.fasta')
parser.add_argument('-c', '--cds', dest = 'c', type = str, help = 'fasta file for the coding region')
parser.add_argument('-s', '--select', dest = 's', type = str, help = 'selected chromosome for analysis')
parser.add_argument('-m', '--model', dest = 'm', type = str, help = 'either a folder of the model files or a model file (keras format) for prediction')
args = parser.parse_args()
genome_fa = pysam.FastaFile(args.f)
vcf_in = VariantFile(args.i)
if args.s:
    params['chr_sele'] = args.s

idf = args.i.split('/')[-1].split('.vcf')[0]
f_log = make_log_file(idf + '_TL_INN_predictions.log', p_params = params, p_vars = vars(args))

t_start = time()

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
pwrite(f_log, '\nChromosome selected in this analysis: ' + params['chr_sele'] + timer(t_start))
bed_dict = {}
with open(args.b, 'r') as f:
    for line in f:
        lst = line.strip().split('\t')
        chromosome = lst[0]
        if chromosome == params['chr_sele']:
            pa = lst[3]
            if pa in bed_dict:
                bed_dict[pa].append(lst)
            else:
                bed_dict[pa] = [lst]
pwrite(f_log, 'Number of poly(A) sites in chromosome ' + params['chr_sele'] + ': ' + str(len(bed_dict)))
if len(bed_dict) == 0:
    stop('No poly(A) sites in the input bed file. Maybe wrong chromosome name?')

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

######------------------------------------------------------------------------------------------------------------
# iterate through all poly(A) sites
pwrite(f_log, '\nStart predicting poly(A) tail length of the variants...' + timer(t_start))
feed_counter = 1
chunk_percent = 10
out_df = pd.DataFrame(columns=['chromosome', 'position', 'ref', 'alt', 'af', 'ac', 'an', 'id', 'type', 'dis2end', 'cpe', 'pas', 'ref_pa', 'ref_tl', 'diff_tl'])
for n_count, pa in enumerate(bed_dict):
    if (len(bed_dict[pa]) > 1):
        if bed_dict[pa][0][5] == '+':
            bed_dict[pa] = sorted(bed_dict[pa], key = lambda x: int(x[1]))
        else:
            bed_dict[pa] = sorted(bed_dict[pa], key = lambda x: int(x[1]), reverse = True)

    # get the reference 3'-UTR sequence
    wt_seq_lst = []
    for i in range(len(bed_dict[pa])):
        if bed_dict[pa][i][5] == '+':
            seq = genome_fa.fetch(bed_dict[pa][i][0], int(bed_dict[pa][i][1]), int(bed_dict[pa][i][2]))
        else:
            seq = reverse_complement(genome_fa.fetch(bed_dict[pa][i][0], int(bed_dict[pa][i][1]), int(bed_dict[pa][i][2])))
        wt_seq_lst.append(seq)
    wt_seq = ''.join(wt_seq_lst)

    # get the CDS sequence if necessary
    if params['flag_cds_track']:
        transcript_id = pa.split('__')[1]
        if transcript_id in cds_dict:
            cds_seq = cds_dict[transcript_id]
        else:
            cds_seq = ''

    # 3'-UTR must be at least this long for analysis
    if len(wt_seq + params['seq_append']) >= params['len_min']: 

        # for each exon region, get the variant info, and construct the entire 3'-UTR sequence associated with each variant
        variant_seq_lst = []
        variant_df = pd.DataFrame(columns=['chromosome', 'position', 'ref', 'alt', 'af', 'ac', 'an', 'id', 'type', 'dis2end', 'cpe', 'pas', 'exon_info', 'ref_pa'])
        for exon_idx in range(len(bed_dict[pa])):
            for rec in vcf_in.fetch(bed_dict[pa][exon_idx][0], int(bed_dict[pa][exon_idx][1]), int(bed_dict[pa][exon_idx][2])):
                exon_info = '_'.join([bed_dict[pa][exon_idx][j] for j in [0,1,2,5]])
                
                # get basic infomation
                v_chr = rec.chrom
                v_pos = rec.pos # 1-based from a VCF file
                v_ref = rec.ref
                if 'AN' in rec.info.keys():
                    v_an = rec.info['AN']
                else:
                    v_an = 'NA'
                
                if rec.id is None:
                    v_id = 'NA'
                else:
                    v_id = rec.id

                if 'variant_type' in rec.info.keys():
                    v_type = rec.info['variant_type']
                else:
                    v_type = 'NA'

                # each allele may have multiple alternative alleles
                for allele_idx in range(len(rec.alts)):

                    # only analyze variants that pass all filters
                    if list(map(str, rec.filter))[allele_idx] == 'PASS':
                        v_alt = rec.alts[allele_idx]
                        if 'AF' in rec.info.keys():
                            v_af = rec.info['AF'][allele_idx]
                        else:
                            v_af = 'NA'
                        if 'AC' in rec.info.keys():
                            v_ac = rec.info['AC'][allele_idx]
                        else:
                            v_ac = 'NA'
                        
                        # get the variant 3'-UTR sequence 
                        # Note that in cases when the ref allele is a long deletion, it can start outside the exon, leads to "v_pos < int(bed_dict[pa][i][1])"
                        if bed_dict[pa][exon_idx][5] == '+':
                            ref_exon_seq = wt_seq_lst[exon_idx]
                            alt_exon_seq = ref_exon_seq[:max((v_pos - int(bed_dict[pa][exon_idx][1]) - 1),0)] + v_alt + ref_exon_seq[(v_pos - int(bed_dict[pa][exon_idx][1]) - 1 + len(v_ref)):]
                        else:
                            ref_exon_seq = reverse_complement(wt_seq_lst[exon_idx])
                            alt_exon_seq = ref_exon_seq[:max((v_pos - int(bed_dict[pa][exon_idx][1]) - 1),0)] + v_alt + ref_exon_seq[(v_pos - int(bed_dict[pa][exon_idx][1]) - 1 + len(v_ref)):]
                            alt_exon_seq = reverse_complement(alt_exon_seq)

                        alt_seq_lst = wt_seq_lst[:]
                        alt_seq_lst[exon_idx] = alt_exon_seq
                        alt_seq = ''.join(alt_seq_lst)
                        variant_seq_lst.append(alt_seq)

                        # calculate the distance between the variant position and the 3'-end
                        # note that this metric is meaningful only when the variant type is SNV
                        if bed_dict[pa][exon_idx][5] == '+':
                            dis2end = max(int(bed_dict[pa][exon_idx][2]) - v_pos - len(v_ref) + 1, 0) + np.sum([len(seq) for n in range(len(alt_seq_lst)) if n > exon_idx])
                        else:
                            dis2end = v_pos - int(bed_dict[pa][exon_idx][1]) - 1 + len(v_ref) - 1 + np.sum([len(seq) for n in range(len(alt_seq_lst)) if n > exon_idx])

                        # check if the varaint has eliminated or created new CPE or PAS motifs
                        n_cpe = len(get_motif_positions(alt_seq, 'TTTTA')) - len(get_motif_positions(wt_seq, 'TTTTA'))
                        n_pas = len(get_motif_positions(alt_seq, 'AWTAAA')) - len(get_motif_positions(wt_seq, 'AWTAAA'))          

                        # add information to the metadata dataframe
                        v_info_lst = [v_chr, v_pos, v_ref, v_alt, v_af, v_ac, v_an, v_id, v_type, dis2end, n_cpe, n_pas, exon_info, bed_dict[pa][exon_idx][3]]
                        variant_df = pd.concat([variant_df, pd.DataFrame([v_info_lst], columns=variant_df.columns)], ignore_index=True)

        # convert sequences to 1-hot encoding and predict tail length with the model
        # the last entry is the wild-type sequence
        if len(variant_seq_lst) > 0:
            X_lst = []
            for seq in (variant_seq_lst + [wt_seq]):
                seq = seq + params['seq_append']
                if params['flag_cds_track']:
                    utr_seq = seq
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
        
            y_lst = []
            for model in models:
                if params['flag_multi_heads']:
                    predictions = model.predict([X_mat, np.repeat(params['model_index'], X_mat.shape[0])], verbose = 0, batch_size = params['batch_size']).ravel()
                else:
                    predictions = model.predict(X_mat, verbose = 0, batch_size = params['batch_size']).ravel()

                y_lst.append(predictions)

            y_avg = np.mean(np.stack(y_lst, axis = 0), axis = 0)
            y_wt = y_avg[-1]
            y_diff = y_avg[:-1] - y_wt
            variant_df['ref_tl'] = y_wt
            variant_df['diff_tl'] = y_diff

            out_df = pd.concat([out_df, variant_df], ignore_index=True)
            pwrite(f_log, 'Finished analyzing this pA site (' + str(len(variant_seq_lst)) + ' variants): ' + bed_dict[pa][i][3])
        
        else:
            pwrite(f_log, 'Skipped. No variants found: ' + bed_dict[pa][i][3])

    else:
        pwrite(f_log, 'Skipped. UTR shorter than ' + str(params['len_min']) + ' nt: ' + bed_dict[pa][i][3])

    if (n_count + 1) >= len(bed_dict) * chunk_percent / 100.0 * feed_counter:
        pwrite(f_log, '\n' + str(feed_counter * chunk_percent) + ' percent poly(A) sites (n=' + str(n_count+1) + ') have been processed...' + timer(t_start) + '\n')
        feed_counter += 1

pwrite(f_log, 'Total number of variants processed: ' + str(out_df.shape[0]))
f_log.close()
out_df.to_csv(params['chr_sele'] + '_INN_predictions.txt', index = False, sep = '\t')
            

