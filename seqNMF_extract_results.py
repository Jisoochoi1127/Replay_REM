#%% Imports
import numpy as np
from seqnmf import seqnmf
#from utils.utils import test_significance
import os
from scipy.stats import skew
from tqdm import tqdm
import pandas as pd
import itertools

import matplotlib.pyplot as plt
plt.style.use('plot_style.mplstyle')
import scipy.io as sio
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
    for shuffle_i in range(params['numShuffles']):
        shuffled_W = np.zeros(W_train.shape)
        for neuron in range(params['numNeurons']):
            shuffled_W[neuron,:,:] = np.roll(W_train[neuron,:,:],shift=np.random.randint(W_train.shape[2]),axis=1)
        
        shuffled_test_H = extract_H(shuffled_W,test_data)
        seq_shuffled_score[shuffle_i,:] = skew(shuffled_test_H,axis=1)

    seq_zscore = np.zeros(params['K']) # zscore for each sequences
    seq_pvalue = np.zeros(params['K']) # zscore for each sequences
    for k in range(params['K']):
        # Z-score
        seq_zscore[k] = (seq_score[k]-np.mean(seq_shuffled_score[:,k]))/np.std(seq_shuffled_score[:,k])

        # p-value
        seq_pvalue[k] = sum(seq_shuffled_score[:,k]>seq_score[k])/params['numShuffles']
    
    return seq_score, seq_shuffled_score, seq_zscore ,seq_pvalue

def extract_seqReplay_score(data_ref, data_pred, params):

    # Extract sequences
    W_ref, H_ref, _, _, _ = seqnmf(
        data_ref.T,
        K=params['K'],
        L=params['L'],
        Lambda=params['Lambda'],
        max_iter=params['maxIters']
        )

    # Extract sequences on predicted data set
    H_pred = extract_H(W_ref, data_pred)
    # output will have d dimensions corresponding to each sequence template
    seqReplay_score = skew(H_pred,axis=1) 

    # Shuffle
    seq_shuffled_score = np.zeros((params['numShuffles'], params['K']))
    for shuffle_i in range(params['numShuffles']):
        shuffled_W = np.zeros(W_ref.shape)
        for neuron in range(params['numNeurons']):
            shuffled_W[neuron,:,:] = np.roll(W_ref[neuron,:,:],shift=np.random.randint(W_ref.shape[2]),axis=1)
        
        shuffled_pred_H = extract_H(shuffled_W,data_pred)
        seq_shuffled_score[shuffle_i,:] = skew(shuffled_pred_H,axis=1)

    seqReplay_zscore = np.zeros(params['K']) # zscore for each sequences
    seqReplay_pvalue = np.zeros(params['K']) # zscore for each sequences
    for k in range(params['K']):
        # Z-score
        seqReplay_zscore[k] = (seqReplay_score[k]-np.mean(seqReplay_shuffled_score[:,k]))/np.std(seqReplay_shuffled_score[:,k])

        # p-value
        seqReplay_pvalue[k] = sum(seqReplay_shuffled_score[:,k]>seqReplay_score[k])/params['numShuffles']

    #TODO add non-parametric pvalue
    return seqReplay_score, seqReplay_shuffled_score, seqReplay_zscore, seqReplay_pvalue

#%% First, measure 'sequenceness' in each individuate session
seqScore_list = []
for condition, mouse, state in tqdm(list(itertools.product(condition_list,
                                                           mouse_list,
                                                           states_list)),
                                                           total=len(condition_list)*len(mouse_list)*len(states_list)):
    # Load data
    data = load_data(mouse, condition, state, params)

    # Extract seq score
    seq_score, seq_shuffled_score, seq_zscore, seq_pvalue = extract_seq_score(data, params)

    seqScore_list.append(
        {
        'mouse':mouse,
        'condition':condition,
        'state':state,
        'max_seq_score': np.max(seq_score), # Take max possible score
        'max_seq_zscore': np.max(seq_zscore)
    })

df=pd.DataFrame(seqScore_list)
df.to_csv(os.path.join(params['path_to_output'],'seqScores.csv'))

#%% Same but look at replay between conditions
seqReplayScore_list = []
for condition, mouse, state_ref, state_pred in tqdm(list(itertools.product(condition_list,
                                                           mouse_list,
                                                           states_list,
                                                           states_list)),
                                                           total=len(condition_list)*len(mouse_list)*len(states_list)*len(states_list)):

    # Load data for both states
    data_ref = load_data(mouse, condition, state_ref, params)
    data_pred = load_data(mouse, condition, state_pred, params)

    # Extract seq score
    seqReplay_score, seqReplay_shuffled_score, seqReplay_zscore, seqReplay_pvalue = extract_seqReplay_score(data_ref, data_pred, params)

    seqReplayScore_list.append(
        {
        'mouse':mouse,
        'condition':condition,
        'state_ref':state_ref,
        'state_pred':state_pred,
        'max_seq_score': np.max(seqReplay_score), # Take max possible score
        'mean_seqReplay_shuffled_score': np.mean(seqReplay_shuffled_score), # Take max possible score
        'max_seqReplay_zscore': np.max(seqReplay_zscore)
    })

df_replay=pd.DataFrame(seqReplayScore_list)
df_replay.to_csv(os.path.join(params['path_to_output'],'seqReplayScores.csv'))

