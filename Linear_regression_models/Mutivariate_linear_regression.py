# v. 20250320
import sys, time, subprocess, argparse, os, shlex, optuna
import numpy as np
import pandas as pd
import tensorflow as tf
from time import time
from datetime import datetime
from tensorflow import keras
from keras import layers
from sklearn.preprocessing import RobustScaler, MaxAbsScaler, PowerTransformer, StandardScaler
from sklearn.model_selection import train_test_split, StratifiedKFold
from scipy import stats
from collections import OrderedDict
from matplotlib import pyplot
from scipy.stats import gaussian_kde, spearmanr, pearsonr

params_global = OrderedDict([
	('flag_input_header', True), # whether the input file with target values has a header
	('input_id_col', 1), # column number in the input file for IDs
	('input_y_col', 2), # column number in the input file for target values
	('input_filter_col', 3), # column number in the input file for filtering rows
	('input_cutoff', 50), # tag cutoff
	('n_cpe_max', 4), # max number of CPEs for labeling each mRNA (if an mRNA has more than this number of CPEs, the 'label' will be set at this number) 
	#('y_mode', 'none'), # how to transform target (tail length) value, can also be changed by input option '-t'
	#('tl_cons', 35), # tail length constant to be subtracted if 'y_mode' is set to 'diff'
	('test_size', 0.1), # fraction of the dataset used as test set
	('n_rep', 5), # number of repeated training with 1 cross-validation group
	#('n_threads', 72), # number of threads for multiprocessing
	#('chunk', 20000), # number of lines for each batch of multiprocessing
	('verbose', 0),
	('n_trials', 1000) # number of trials for optimization
])

# the following are parameters used for optimization

model_params = OrderedDict([
	('epochs', 100),
	('batch_size', 50), 
	('learning_rate', 0.0025169),
	('l1_reg', 0.000166),
	('l2_reg', 0.0001097), 
	('loss', 'mse'),
	('metrics', 'mae'),
	('patience', 10)  # number of epochs with no improvement before early stopping
])

param_search_dict = OrderedDict([
	('learning_rate', {'min': 0.005, 'max': 0.5}),# {'min': 0.00001, 'max': 0.1}),
	('l1_reg', {'min': 0.000001, 'max': 0.01}),
	('l2_reg', {'min': 0.000001, 'max': 0.01}),
	('batch_size', [1000, 2000, 3000, 7000]), #[50, 100, 200, 500, 1000]),
	('loss', ['mse']), #['mse', 'mae']),
	('epochs', [30, 50, 100])
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

def fasta_id_parser(line, sep = '\t', pos = 0):
	sele_id = line.rstrip().lstrip('>').split(sep)[pos]
	return(sele_id)

def fasta_to_dict(file, len_min = 0, string_filter = None):
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
						current_id = fasta_id_parser(line, sep = '__', pos = 0)
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

def model_linear_regression(x_train, y_train, x_val, y_val, params, best_model = False, early_stop = False):
	'''
	Train a linear regression model with the specified parameters.

	Args:
		x_train, y_train: Training data and target values.
		x_val, y_val: Validation data and target values.
		params: Dictionary of hyperparameters.
		best_model: If True, loads the best saved model after training.
		early_stop: If True, uses early stopping.

	Returns:
		Tuple (training history, trained model, updated params).
	'''

	# Ensure global variables are initialized
	global params_global, model_params, f_log, args, idf
	
	# update missing parameters with defaults
	for key in model_params:
		if key not in params:
			params.update({key: model_params[key]})

	# define the model
	model = keras.Sequential([
		layers.Input(shape=(x_train.shape[1],)),  
		layers.Dense(1, kernel_regularizer = keras.regularizers.l1_l2(l1=params['l1_reg'], l2=params['l2_reg']))
	])
	
	model.compile(
		optimizer = tf.optimizers.Adam(learning_rate = params['learning_rate']), 
		loss = params['loss'], 
		metrics = [params['metrics']]
	)

	# Define callbacks
	f_best_model = os.path.join(args.o, f'{idf}_callback_best_model.keras')
	callbacks = [keras.callbacks.ModelCheckpoint(f_best_model, monitor = 'val_loss', mode = 'min', save_best_only = True)]
	if early_stop:
		callbacks.append(keras.callbacks.EarlyStopping(monitor = 'val_loss', mode = 'min', patience = params['patience']))

	history = model.fit(
		x = x_train, y = y_train, 
		batch_size = params['batch_size'], 
		epochs = params['epochs'], 
		verbose = params_global['verbose'], 
		shuffle = True, 
		callbacks = callbacks, 
		validation_data=(x_val, y_val)
	)

	if best_model:
		pwrite(f_log, 'Return the best model...\n')
		try:
			model = keras.models.load_model(f_best_model)
			#model = keras.models.load_model(f_best_model, custom_objects={'custom_metric': inn_models.custom_metric})
		except Exception as e:
			pwrite(f_log, f"Error loading best model: {e}\n")
			pwrite(f_log, 'Return the final model...\n')
	else:
		pwrite(f_log, 'Return the final model...\n')

	return(history, model, params)


def optuna_objective(trial, x_train, y_train, x_val, y_val, params, opti_log = None):
	'''
	Optuna objective function for hyperparameter optimization.

	Args:
		trial: Optuna trial object.
		x_train, y_train: Training data and target values.
		x_val, y_val: Validation data aand target values.
		params: Dictionary of hyperparameters.
		opti_log: Path to the log file for optimization results.

	Returns:
		Negative R-squared value for Optuna optimization.
	'''

	# Ensure global variables are initialized
	global params_global, args, idf, best_neg_r2, f_log
	
	# Sample hyperparameters
	params_opti = OrderedDict()
	for param in params:
		if param in ['l1_reg', 'l2_reg', 'learning_rate']:
			if params[param]['min'] < params[param]['max']:
				params_opti[param] = trial.suggest_float(param, params[param]['min'], params[param]['max'], log = True)
		else:
			if isinstance(params[param], list) and len(params[param]) > 1:
				params_opti[param] = trial.suggest_categorical(param, params[param])
			else:
				params_opti[param] = params[param][0] if isinstance(params[param], list) else params[param]

	# Train model using these hyperparameters
	history, model, out_params = model_linear_regression(x_train, y_train, x_val, y_val, params_opti, best_model = True, early_stop = True)

	# Predict using the trained model
	y_pred = model.predict(x_val).ravel()

	# Compute Negative R-squared for minimization
	neg_r2 = -stats.linregress(y_val,y_pred).rvalue ** 2	

	# Save trial results
	if opti_log:
		if not os.path.exists(opti_log):	
			with open(opti_log, 'w') as f:
				f.write('Trial number\t' + '\t'.join(list(trial.params.keys())) + '\tepochs (actual/max)\tneg_r2\n')
				f.write(f'{trial.number + 1}\t' + '\t'.join([str(trial.params[param]) for param in trial.params]) + 
					f"\t{len(history.history['loss'])}/{out_params['epochs']}\t{neg_r2}\n")
		else:
			with open(opti_log, 'a') as f:	
				f.write(f'{trial.number + 1}\t' + '\t'.join([str(trial.params[param]) for param in trial.params]) + 
					f"\t{len(history.history['loss'])}/{out_params['epochs']}\t{neg_r2}\n")

	if neg_r2 < best_neg_r2:
		best_neg_r2 = neg_r2

		# Generate scatter plot
		scatter_plot(y_val, y_pred, fn = os.path.join(args.o, 'Optimization_best_model_predicted_vs_measured_scatter_plot.png'))
		
	pwrite(f_log, f"\nCurrent metric for trial number {trial.number + 1} is {neg_r2:.3f}. Best metric so far is {best_neg_r2:.3f}.\n")
	return(neg_r2)

def save_opti_results(study, fn = 'Optimization_reuslts'):
	study.trials_dataframe().to_csv(f'{fn}_results.csv')

	# make a few plots
	# optimization history
	plt = optuna.visualization.matplotlib.plot_optimization_history(study).figure
	plt.savefig(f'{fn}_optimization_history.png', dpi = 300, bbox_inches='tight')

	# param slices
	plt = optuna.visualization.matplotlib.plot_slice(study)[0].figure
	plt.savefig(f'{fn}_search_slice.png', dpi = 300, bbox_inches='tight')

	# param importances
	plt = optuna.visualization.matplotlib.plot_param_importances(study).figure
	plt.savefig(f'{fn}_param_importances.png', dpi = 300, bbox_inches='tight')

	# timeline
	plt = optuna.visualization.matplotlib.plot_timeline(study).figure
	plt.savefig(f'{fn}_timeline.png', dpi = 300, bbox_inches='tight')

	# contour plots
	n_params = len(list(study.best_trial.params.keys()))
	fig_size = max(3, n_params * 3)
	plt = optuna.visualization.matplotlib.plot_contour(study)[0,0].figure
	plt.set_size_inches(fig_size, fig_size) 
	plt.savefig(f'{fn}_contour_map.png', dpi = 300, bbox_inches='tight')


def scatter_plot(x, y, fn):

	# Calculate point density
	xy = np.vstack([x, y])
	z = gaussian_kde(xy)(xy)

	# Calculate Spearman and Pearson correlation coefficients
	spearman_corr, _ = spearmanr(x, y)
	pearson_corr, _ = pearsonr(x, y)

	# Calculate the number of points
	n = len(x)

	# Create scatter plot
	pyplot.figure(figsize=(8, 6))
	scatter = pyplot.scatter(x, y, c=z, cmap='viridis', s=20)

	# Add labels
	pyplot.xlabel('Measured data')
	pyplot.ylabel('Predicted data')

	# Add color bar
	pyplot.colorbar(scatter, label='Density')

	# Ensure x and y axes have the same scale
	pyplot.axis('equal')

	# Plot dashed line at x = y
	min_val = min(min(x), min(y))
	max_val = max(max(x), max(y))
	pyplot.plot([min_val, max_val], [min_val, max_val], '--', color='red', linewidth=1)

	# Add correlation and sample size text in the bottom right corner
	pyplot.text(0.95, 0.05, f'n={n}\nPearson: {pearson_corr:.2f}\nSpearman: {spearman_corr:.2f}',
		horizontalalignment='right', verticalalignment='bottom',
		transform=pyplot.gca().transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

	# Save the plot to a PNG file
	pyplot.savefig(fn, dpi=300, bbox_inches='tight')

######------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest = 'i', type = str, required = True, help = 'input file name')
parser.add_argument('-u', '--utr', dest = 'u', type = str, required = True, help = 'fasta file')
parser.add_argument('-f', '--feature', dest = 'f', type = str, help = 'info file for files with features')
parser.add_argument('-o', '--output', dest = 'o', type = str, default = './', help = 'output folder')
parser.add_argument('-m', '--mode', dest = 'm', type = str, default = 'test', help = 'mode to run this script, can be "cv", "opti", "test", or "predict"')
parser.add_argument('-op', '--hyperparameters', dest = 'op', type = int, help = 'Number of hyperparameter combinations')

#parser.add_argument('-t', '--transform', dest = 't', type = str, default = 'none', help = 'how to transform Y, can be "none", "diff", log", "sqrt", "boxcox", "yeo-johnson"')
args = parser.parse_args()

t_start = time() # timer start

os.makedirs(args.o, exist_ok=True)

if args.op:
	params_global['n_trials'] = args.op

# construct out file prefix based on different input files
idf = args.i.split('/')[-1].split('.')[0] 
f_log = make_log_file(os.path.join(args.o, idf + '_log.txt'), p_params = params_global, p_vars = vars(args))

pwrite(f_log, 'Linear regression with Adam optimizer:' + timer(t_start))

# make a dictionary of confident UTRs
pwrite(f_log, '\n' + ('').join(['-']*100))
pwrite(f_log, 'Make a dictionary of UTR sequences with confident poly(A) sites...' + timer(t_start))
utr_dict = fasta_to_dict(args.u, len_min = 0, string_filter = 'uncertain')
pwrite(f_log, 'Number of entries in the fasta file after filtering: ' + str(len(utr_dict)))

# get y target values
pwrite(f_log, '\n' + ('').join(['-']*100))
pwrite(f_log, 'Obtain target values...' + timer(t_start))
id_lst = []
y_lst = []
with open(args.i, 'r') as f:
	if params_global['flag_input_header']:
		f.readline()
	for i, line in enumerate(f):
		lst = line.strip().split('\t')
		if (params_global['input_cutoff'] is None or float(lst[params_global['input_filter_col']-1]) >= params_global['input_cutoff']) and \
		lst[params_global['input_y_col']-1] != 'NA' and lst[params_global['input_id_col']-1] in utr_dict:
			id_lst.append(lst[params_global['input_id_col']-1])
			y_lst.append(float(lst[params_global['input_y_col']-1]))
	y_df = pd.DataFrame({'id':id_lst, 'y':y_lst})

pwrite(f_log, 'Number of entries in the input file: ' + str(i+1))
pwrite(f_log, 'Number of entries passing the cutoff: ' + str(len(y_df)))


# get pre-calcuated features
pwrite(f_log, '\n' + ('').join(['-']*100))
pwrite(f_log, 'Obtain feature values from files...' + timer(t_start))
df_ff = pd.read_csv(args.f, delimiter='\t')
df_feature = pd.DataFrame({'id':id_lst})
for i in range(df_ff.shape[0]):
	pwrite(f_log, str(df_feature.shape))
	feature_type = df_ff['type'].iloc[i]
	feature_file = df_ff['file'].iloc[i]
	pwrite(f_log, 'File ' + str(i+1) + ' for ' + str(feature_type) + ':')
	pwrite(f_log, feature_file)
	df_temp = pd.read_csv(feature_file, delimiter='\t')
	df_feature = pd.merge(df_feature, df_temp, on = 'id', how = 'inner')
pwrite(f_log, 'Number of entries: ' + str(df_feature.shape[0]))
pwrite(f_log, 'Number of features: ' + str(df_feature.shape[1] - 1))

# take the common entries of the target and the feature dataframes and make sure they have the same order
y_df = pd.merge(y_df, df_feature[['id', 'num_TTTTA']], on = 'id', how = 'inner')
y_df['label'] = np.where(y_df['num_TTTTA'] <= params_global['n_cpe_max'], y_df['num_TTTTA'], 4) # add CPE label to each entry 
y_df = y_df.drop('num_TTTTA', axis = 1)
df_feature = pd.merge(y_df[['id']], df_feature, on = 'id', how = 'inner')
x_df = df_feature.drop('id', axis = 1)

pwrite(f_log, '\nData dimension, X: ' + str(x_df.shape) + ' Y: ' + str(y_df.shape))

###---------------------------------------------------------------------------------------------------------
# Set available GPUs
pwrite(f_log, '\n' + ('').join(['-']*100))
gpus = tf.config.list_physical_devices('GPU')
if gpus:
	pwrite(f_log, f'Detected {len(gpus)} GPU cores.')
	n_jobs = len(gpus)
else:
	pwrite(f_log, f'No GPU core is detected. Use CPU for this session. ')
	n_jobs = 1

###---------------------------------------------------------------------------------------------------------
# regression

# for optimization of hyperparameters
if args.m.startswith('op'):
	pwrite(f_log, '\n' + ('').join(['-']*100))
	pwrite(f_log, 'Performing hyperparameter search...'+ timer(t_start))
	pwrite(f_log, 'Ratio between testing and training sets: ' + str(params_global['test_size']))

	# split data into training and testing groups
	x_train, x_test, y_train, y_test= train_test_split(x_df, y_df, test_size = params_global['test_size'], random_state = 57, shuffle = True, stratify = y_df['label'])
	y_train = y_train['y'].to_numpy()
	y_test = y_test['y'].to_numpy()
	
	# scale data
	x_scaler = RobustScaler()
	x_scaler.fit(x_train)
	x_train_scaled = x_scaler.transform(x_train)
	x_test_scaled = x_scaler.transform(x_test)

	#y_scaler = PowerTransformer()
	#y_scaler.fit(y_train.reshape(-1,1))
	#y_train_trans = y_scaler.transform(y_train.reshape(-1,1)).ravel()
	#y_test_trans = y_scaler.transform(y_test.reshape(-1,1)).ravel()
	#y_train_trans = y_train[:]
	#y_test_trans = y_test[:]

	# optimize with objective function
	best_neg_r2 = 0
	f_opti = os.path.join(args.o, f'{idf}_hyperparameter_search_log.txt')
	study_storage = 'sqlite:///' + os.path.join(args.o, f'{idf}_hyperparameter_search_optuna_study.db')
	
	try:
		optuna.delete_study(study_name = 'Linear_regression', storage = study_storage)
		pwrite(f_log, '\nAn optuna study with the same name exists. This old one will be overwritten by the new study.')
	except KeyError:
		pwrite(f_log, '\nA new optuna study is created for optimization.')

	study = optuna.create_study(
		study_name = 'Linear_regression',
		direction = 'minimize', 
		sampler = optuna.samplers.TPESampler(),
		storage = study_storage
		)
	
	study.optimize(
		lambda trial: optuna_objective(trial, x_train = x_train_scaled, y_train = y_train, x_val = x_test_scaled, y_val = y_test, params = param_search_dict, opti_log = f_opti), 
		n_trials = params_global['n_trials'], n_jobs = 1, gc_after_trial = True
		)

	# save optimization results
	save_opti_results(study = study, fn = os.path.join(args.o, f'{idf}_hyperparameter_search'))

# for training and testing with defined hyperparameters
elif args.m == 'cv' or args.m == 'test':
	pwrite(f_log, '\n' + ('').join(['-']*100))
	pwrite(f_log, 'Use a linear regression model for training and testing in a CV fold of ' + str(int(1/params_global['test_size'])) + ' with the following hyperparameters...' + timer(t_start))
	
	# split data into k-fold, stratified by "label"
	sfk_split = StratifiedKFold(n_splits = int(1/params_global['test_size']), shuffle = True, random_state = 57).split(x_df, y_df['label'])

	out_pred = os.path.join(args.o, idf + '_train_test_CV_prediction_results.txt')
	coef_df = pd.DataFrame({'feature': x_df.columns.to_list() + ['intercept']})
	with open(out_pred, 'w') as f_res, open(os.path.join(args.o, idf + '_train_test_CV_stats.txt'), 'w') as f_stat:
		f_stat.write('cv_group\trep\tr2\tr_pearson\tr_spearman\tloss\tmetric\tepochs\n')
		f_res.write('group\tid\ty\ty_pred_avg\ty_pred_best\tlabel\n')
		for i, (train_val_idx, test_idx) in enumerate(sfk_split):
			r2_lst = []
			
			x_train_val = x_df.loc[train_val_idx].to_numpy()
			x_test = x_df.loc[test_idx].to_numpy()
			y_train_val = y_df.loc[train_val_idx]['y'].to_numpy()
			y_test = y_df.loc[test_idx]['y'].to_numpy()

			# scale data and transform to tensors
			x_scaler = RobustScaler()
			x_scaler.fit(x_train_val)
			x_train_val_scaled = x_scaler.transform(x_train_val)
			x_test_scaled = x_scaler.transform(x_test)

			for j in range(params_global['n_rep']):
				# train model
				pwrite(f_log, '\nTraining with CV group ' + str(i+1) + ' out of ' + str(int(1/params_global['test_size'])) + ', replicate ' + str(j+1) + ' out of ' + str(params_global['n_rep']) + timer(t_start))
				
				# split data into train and validation sets
				x_train, x_val, y_train, y_val = train_test_split(x_train_val_scaled, y_train_val, test_size = params_global['test_size'], random_state = i*10+j, shuffle = True, stratify = y_df.loc[train_val_idx]['label'].to_numpy())
				
				history, model, out_params = model_linear_regression(x_train, y_train, x_val, y_val, model_params, best_model = True, early_stop = True)
				
				# predict the test set and evaluate the model:
				scores = model.evaluate(x_test_scaled, y_test, batch_size = model_params['batch_size'], verbose = params_global['verbose'])
				y_test_pred = model.predict(x_test_scaled, batch_size = model_params['batch_size'], verbose = params_global['verbose']).ravel()
				
				r_pearson = stats.pearsonr(y_test, y_test_pred).statistic
				r_spearmam = stats.spearmanr(y_test, y_test_pred).statistic
				r2_lst.append(r_pearson ** 2)
				pwrite(f_log, f'R squared: {r_pearson ** 2: .2f} (Rp = {r_pearson: .3f}, Rs = {r_spearmam: .3f})\n')

				# output stats
				f_stat.write(
					f'{i+1}\t{j+1}\t{r_pearson ** 2: .2f}\t{r_pearson: .3f}\t{r_spearmam: .3f}\t' + 
					'\t'.join(list(map(str, scores))) + 
					f'\t{len(history.history["loss"])}/{model_params["epochs"]}\n')
				
				if args.m != 'cv': # only do this once if not in 'cv' mode
					break

				# get the weights
				weights, bias = model.layers[0].get_weights()
				coef_df[f"cv_{i+1}_rep_{j+1}"] = np.append(weights.flatten(), bias).ravel()

				# combine predicted values from each rep for each cv split
				if j == 0:
					y_pred_cv = np.array(y_test_pred)
				else:
					y_pred_cv = np.vstack((y_pred_cv, np.array(y_test_pred)))

			if args.m != 'cv': # only do this once if not in 'cv' mode
				break
			
			# average predicted values among all reps for each cv split
			y_pred_avg = np.mean(y_pred_cv, axis = 0)
			y_pred_best = y_pred_cv[r2_lst.index(max(r2_lst)), :]
			
			for k in range(len(y_test)):
				f_res.write(f"{i+1}\t{y_df.loc[test_idx]['id'].to_numpy()[k]}" + 
					f"\t{y_test[k]}\t{y_pred_avg[k]}\t{y_pred_best[k]}\t{y_df.loc[test_idx]['label'].to_numpy()[k]}\n")

	# output the coefficient 
	out_coef = os.path.join(args.o, 'Linear_Regression_coefficient_' + idf + '_by_CV.csv')
	coef_df.to_csv(out_coef, index=False)
		
	# Make two plosts (require an R script "MPRA_linear_regression_plots.R")
	# 1. A scatter plot for comparing measured and predicted y values
	# 2. A heatmap for coefficients of top position kmers
	if os.path.exists('Linear_regression_plots.R'):
		pwrite(f_log, '\nMaking plots with an R script...')
		command = 'Rscript Linear_regression_plots.R ' + out_pred + ' ' + out_coef
		subprocess.Popen(shlex.split(command)).communicate()
else: # args.m == 'predict'
	pass

# Make a plot for comparing measured and predicted values 
if 'out_pred' in globals() and os.path.exists(out_pred):
	out_df = pd.read_csv(out_pred, sep = '\t')
	scatter_plot(out_df['y'], out_df['y_pred_avg'], fn = os.path.join(args.o, 'idf' + '_pred_avg_final_scatter_plot.png'))
	scatter_plot(out_df['y'], out_df['y_pred_best'], fn = os.path.join(args.o, 'idf' + '_pred_best_final_scatter_plot.png'))

pwrite(f_log, 'Finished...' + timer(t_start))
f_log.close()


