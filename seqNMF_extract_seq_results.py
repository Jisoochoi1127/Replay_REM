#%% Imports
import numpy as np
from seqnmf import seqnmf
import os
from scipy.stats import skew
from tqdm import tqdm
import pandas as pd
import itertools
import scipy.io as sio
from scipy.signal import find_peaks
import yaml
import h5py
from pycaan.functions.signal_processing import binarize_ca_traces

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

#%% Define conditions, dataset
states_list = ['REMpre', 'wake', 'REMpost']
condition_list = ['LTD1','LTD5','HATD5']
mouse_list = ['pv1060', 'pv1254', 'pv1069']

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
    W_train, _, _, _, _ = seqnmf(
        train_data.T,
        K=params['K'],
        L=params['L'],
        Lambda=params['Lambda'],
        max_iter=params['maxIters']
        )

    # Extract sequences on test set
    H_test = extract_H(W_train, test_data)
    # output will have d dimensions corresponding to each sequence template
    seq_scores = skew(H_test,axis=1) 

    # Shuffle
    seq_shuffled_score = np.zeros((params['numShuffles'], params['K']))
    shuffled_test_H = np.zeros((params['numShuffles'], params['K'],len(test_data)))

    for shuffle_i in range(params['numShuffles']):
        shuffled_W = np.zeros(W_train.shape)
        for neuron in range(params['numNeurons']):
            shuffled_W[neuron,:,:] = np.roll(W_train[neuron,:,:],shift=np.random.randint(W_train.shape[2]),axis=1)
        
        temp = extract_H(shuffled_W,test_data)
        shuffled_test_H[shuffle_i,:,:] = temp
        seq_shuffled_score[shuffle_i,:] = skew(temp,axis=1)

    seq_pvalues = np.zeros(params['K']) # zscore for each sequences
    H_test_pvalue = np.zeros((params['K'],len(test_data)))
    zscored_H = np.zeros((params['K'],len(test_data)))
    H_test_confidence = np.zeros((params['K'],len(test_data)))
    seq_locs = {}

    for k in range(params['K']):
        # z-score actual H to extract signal-to-noise:
        zscored_H[k,:] = (H_test[k,:]-np.mean(H_test[k,:]))/np.std(H_test[k,:])
        for i in range(len(test_data)):
            H_test_pvalue[k,i] = sum(shuffled_test_H[:,k,i]>H_test[k,i])/params['numShuffles']
        
        H_test_confidence[k,:] = 1-H_test_pvalue[k,:]
        H_test_confidence[k,zscored_H[k,:]<2] = 0

        peaks, _ = find_peaks(H_test_confidence[k], 
                        height=(1,None), # Only include peaks for confidence==1
                        distance=params['L'])
        seq_locs.update({
            k:peaks
        })

        seq_pvalues[k] = sum(seq_shuffled_score[:,k]>seq_scores[k])/params['numShuffles']
    
    return seq_scores, seq_pvalues, seq_locs

#%% First, measure 'sequenceness' in each individuate session
for condition, mouse, state in tqdm(list(itertools.product(condition_list,
                                                           mouse_list,
                                                           states_list)),
                                                           total=len(condition_list)*len(mouse_list)*len(states_list)):
    
    if not os.path.exists(os.path.join(params['path_to_output'],f'seqResults_{condition}_{mouse}_{state}.h5')):
        # Load data
        data = load_data(mouse, condition, state, params)

        # Extract seq score
        seq_scores, seq_pvalues, seq_locs = extract_seq_score(data, params)

        with h5py.File(os.path.join(params['path_to_output'],f'seqResults_{condition}_{mouse}_{state}.h5'),'w') as f:
            f.create_dataset('mouse', data=mouse)
            f.create_dataset('condition', data=condition)
            f.create_dataset('state', data=state)
            f.create_dataset('S1_score', data=seq_scores[0])
            f.create_dataset('S2_score', data=seq_scores[1])
            f.create_dataset('S1_pvalue', data=seq_pvalues[0])
            f.create_dataset('S2_pvalue', data=seq_pvalues[1])
            f.create_dataset('S1_numSeqs', data=len(seq_locs[0]))
            f.create_dataset('S2_numSeqs', data=len(seq_locs[1]))