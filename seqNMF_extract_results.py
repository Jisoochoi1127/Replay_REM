#%% Imports
import numpy as np
from seqnmf import seqnmf
#from utils.utils import test_significance
import os
from scipy.signal import find_peaks
from scipy.stats import skew
from tqdm import tqdm
import pandas as pd

import matplotlib.pyplot as plt
plt.style.use('plot_style.mplstyle')
import scipy.io as sio
import yaml
import h5py
from pycaan.functions.dataloaders import load_data
from pycaan.functions.signal_processing import binarize_ca_traces, preprocess_data

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

#%% Custom functions
def open_file(path, filename):
    data={}
    try:
        f = h5py.File(os.path.join(path,filename)+'.mat', 'r')
        data.update(
        {
        'rawData':f[filename][()].T
        }
        )
    except:
        f = sio.loadmat(path+'/ms.mat')
        data.update(
        {
        'rawData':f['ms']['RawTraces'][0][0] 
        }
        )
    return data

def load_data(mouse, condition, state, params):
    path=os.path.join(params['path_to_dataset'], mouse, condition)
    if 'pre' in state: 
        filename='all_binary_pre_REM'
        data = open_file(path, filename)

    elif 'post' in state:
        filename='all_binary_post_REM'
        data = open_file(path, filename)

    else: # Then, must be wake data in old mat format
        filename='ms'
        data = open_file(path, filename)

    data['binaryData'], _ = binarize_ca_traces(data['rawData'], 2, params['samplingFrequency'])
    data = data['binaryData'][:,0:params['numNeurons']]
    
    return data

# Custom function to extract H given X and W
def extract_H(W,data):
    H = np.zeros((W.shape[1], data.shape[0]),dtype='float')
    for l in range(W.shape[2]):
        X_shifted = np.roll(data,shift=1-l,axis=0) # Convolution through dataset
        H = H + (X_shifted @ W[:, :, l]).T # Matrix multiplication
    return H

def extract_seq_score(data, params):
    # Define train/test datasets
    trainingFrames = np.zeros(len(data), dtype=bool)
    trainingFrames[0:int(params['train_test_ratio']*len(data))] = True
    testingFrames = ~trainingFrames

    # Split dataset in two
    train_data = data[trainingFrames]
    test_data = data[testingFrames]

    # Extract sequences
    W_train, H_train, _, _, _ = seqnmf(
        train_data.T,
        K=params['K'],
        L=params['L'],
        Lambda=params['Lambda'],
        max_iter=params['maxIters']
        )

    # Extract sequences on test set
    H_test = extract_H(W_train, test_data)
    # output will have d dimensions corresponding to each sequence template
    seq_score = skew(H_test,axis=1) 

    # Shuffle
    seq_shuffled_score = np.zeros((params['numShuffles'], params['K']))
    for shuffle_i in tqdm(range(params['numShuffles'])):
        shuffled_W = np.zeros(W_train.shape)
        for neuron in range(params['numNeurons']):
            shuffled_W[neuron,:,:] = np.roll(W_train[neuron,:,:],shift=np.random.randint(W_train.shape[2]),axis=1)
        
        shuffled_test_H = extract_H(shuffled_W,test_data)
        seq_shuffled_score[shuffle_i,:] = skew(shuffled_H_wake,axis=1)

    seq_zscore = np.zeros(params['K']) # zscore for each sequences
    for k in range(params['K']):
        seq_zscore[k] = (seq_score[k]-np.mean(seq_shuffled_score[:,k]))/np.std(seq_shuffled_score[:,k])

    #TODO add non-parametric pvalue
    return seq_score, seq_shuffled_score, seq_zscore

#%% Define conditions, dataset
seqScore_dict = {}
states_list = ['REMpre', 'wake', 'REMpost']
condition_list = ['LTD1','LTD5','HATD5']
mouse_list = ['pv1060', 'pv1254', 'pv1069']

#%% First, measure 'sequenceness' in each individuate session
for condition in condition_list:
    for mouse in mouse_list:
        for state in states_list:
            print(f'Processing {condition}, {mouse}, {state}...')
            # Load data
            data = load_data(mouse, condition, state, params)

            # Extract seq score
            seq_score, seq_shuffled_score, seq_zscore = extract_seq_score(data, params)

            seqScore_dict.update(
                {
                'mouse':mouse,
                'condition':condition,
                'state':state,
                'max_seq_score': np.max(seq_score), # Take max possible score
                'max_seq_zscore': np.max(seq_zscore)
            }
            )

#%% Save results


#%% Imagine matrix w/ 3 d.o.f.: pre, wake, post REM
# For each mouse/condition what is the seqReplay score across all states?
# How deos it compare to shuffled data?



#%%
for condition in condition_list:
    for ref_state in states_list:
        for target_state in states_list:
            for mouse in mouse_list:
                # Load data
                # Preprocess and binarize
                seqReplay_score, seqReplay_shuffled_score, seqReplay_zscore = extract_seqReplay_score()

#%% Save results
pd.to_h5()
