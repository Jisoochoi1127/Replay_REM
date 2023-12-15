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

def extract_seqReplay_score(data_ref, data_pred, params):

    # Extract sequences
    W_ref, _, _, _, _ = seqnmf(
        data_ref.T,
        K=params['K'],
        L=params['L'],
        Lambda=params['Lambda'],
        max_iter=params['maxIters']
        )

    # Extract sequences on predicted data set
    H_pred = extract_H(W_ref, data_pred)
    # output will have d dimensions corresponding to each sequence template
    seqReplay_scores = skew(H_pred,axis=1) 

    # Shuffle
    seqReplay_shuffled_score = np.zeros((params['numShuffles'], params['K']))
    shuffled_pred_H = np.zeros((params['numShuffles'], params['K'],len(data_pred)))

    for shuffle_i in range(params['numShuffles']):
        shuffled_W = np.zeros(W_ref.shape)
        for neuron in range(params['numNeurons']):
            shuffled_W[neuron,:,:] = np.roll(W_ref[neuron,:,:],shift=np.random.randint(W_ref.shape[2]),axis=1)
        
        
        temp = extract_H(shuffled_W,data_pred)
        shuffled_pred_H[shuffle_i,:,:] = temp
        seqReplay_shuffled_score[shuffle_i,:] = skew(temp,axis=1)

    seqReplay_pvalues = np.zeros(params['K']) # zscore for each sequences
    H_pred_pvalue = np.zeros((params['K'],len(data_pred)))
    zscored_H = np.zeros((params['K'],len(data_pred)))
    H_pred_confidence = np.zeros((params['K'],len(data_pred)))
    seqReplay_locs = {}

    for k in range(params['K']):
        # z-score actual H to extract signal-to-noise:
        zscored_H[k,:] = (H_pred[k,:]-np.mean(H_pred[k,:]))/np.std(H_pred[k,:])
        for i in range(len(data_pred)):
            H_pred_pvalue[k,i] = sum(shuffled_pred_H[:,k,i]>H_pred[k,i])/params['numShuffles']
        
        H_pred_confidence[k,:] = 1-H_pred_pvalue[k,:]
        H_pred_confidence[k,zscored_H[k,:]<2] = 0

        peaks, _ = find_peaks(H_pred_confidence[k], 
                        height=(1,None), # Only include peaks for confidence==1
                        distance=params['L'])
        seqReplay_locs.update({
            k:peaks
        })

        seqReplay_pvalues[k] = sum(seqReplay_shuffled_score[:,k]>seqReplay_scores[k])/params['numShuffles']
    
    return seqReplay_scores, seqReplay_pvalues, seqReplay_locs

#%% Same but look at replay between conditions
seqReplayScore_list = []
for condition, mouse, state_ref, state_pred in tqdm(list(itertools.product(condition_list,
                                                           mouse_list,
                                                           states_list,
                                                           states_list)),
                                                           total=len(condition_list)*len(mouse_list)*len(states_list)*len(states_list)):
    if not os.path.exists(os.path.join(params['path_to_output'],f'seqReplayResults_{condition}_{mouse}_{state_ref}_{state_pred}.h5')):
        # Load data for both states
        data_ref = load_data(mouse, condition, state_ref, params)
        data_pred = load_data(mouse, condition, state_pred, params)

        # Extract seq score
        seqReplay_scores, seqReplay_pvalues, seqReplay_locs = extract_seqReplay_score(data_ref, data_pred, params)

        with h5py.File(os.path.join(params['path_to_output'],f'seqReplayResults_{condition}_{mouse}_{state_ref}_{state_pred}.h5'),'w') as f:
            f.create_dataset('mouse', data=mouse)
            f.create_dataset('condition', data=condition)
            f.create_dataset('state_ref', data=state_ref)
            f.create_dataset('state_pred', data=state_pred)
            f.create_dataset('S1_score', data=seqReplay_scores[0])
            f.create_dataset('S2_score', data=seqReplay_scores[1])
            f.create_dataset('S1_pvalue', data=seqReplay_pvalues[0])
            f.create_dataset('S2_pvalue', data=seqReplay_pvalues[1])
            f.create_dataset('S1_numSeqs', data=len(seqReplay_locs[0]))
            f.create_dataset('S2_numSeqs', data=len(seqReplay_locs[1]))
