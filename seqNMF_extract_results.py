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

#%%
results_dict = {}
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

def extract_seq_score(data)
    trainingFrames = np.zeros(len(full_REMpre_data), dtype=bool)
    trainingFrames[0:int(params['train_test_ratio']*len(full_REMpre_data))] = True
    testingFrames = ~REMpre_trainingFrames

    REMpre_train_data = full_REMpre_data[REMpre_trainingFrames]
    REMpre_test_data = full_REMpre_data[REMpre_testingFrames]


    W_REMpre_train, H_REMpre_train, _, _, _ = seqnmf(
        REMpre_train_data.T,
        K=K,
        L=L,
        Lambda=Lambda,
        max_iter=maxIters
        )


    return seq_score, seq_shuffled_score, seq_zscore

#%% First, measure 'sequenceness' in each individuate session
for condition in condition_list:
    for state in states_list:
        for mouse in mouse_list:
            # Load data
            data = load_data(data, params)

            # Extract seq score
            seq_score, seq_shuffled_score, seq_zscore = extract_seq_score(data, params)

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



#%% Custom function to extract H given X and W
def extract_H(W,data):
    H = np.zeros((W.shape[1], data.shape[0]),dtype='float')
    for l in range(W.shape[2]):
        X_shifted = np.roll(data,shift=1-l,axis=0)
        H = H + (X_shifted @ W[:, :, l]).T #Matrix multiplication
    return H

#%% extract H on test data
H_REMpost_from_pre = extract_H(W_REMpre_train, full_REMpost_data)

#%% Plot resulting sequences
plt.figure(figsize=(7,3))
plt.subplot(232)
plt.imshow(REMpre_train_data[:,sorting_index].T,aspect='auto',cmap='magma',vmax=.5)
plt.axis('off')
plt.title('Pre-task REM')

plt.subplot(233)
plt.imshow(full_REMpost_data[:,sorting_index].T,aspect='auto',cmap='magma',vmax=.5)
plt.axis('off')
plt.title('Post-task REM')

plt.subplot(4,3,8)
plt.plot(H_REMpre_train[0],label='S1',color='C0')
plt.plot(H_REMpre_train[1],label='S2',color='C6')
plt.xlim(0,len(H_REMpre_train[0]))
plt.xticks([])
plt.yticks([])

plt.legend(bbox_to_anchor=(0, 1), loc='upper right', borderaxespad=0)

plt.subplot(4,3,9)
plt.plot(H_REMpost_from_pre[0],label='S1',color='C0')
plt.plot(H_REMpost_from_pre[1],label='S2',color='C6')
plt.xlim(0,len(H_REMpost_from_pre[0]))
plt.xticks([])
plt.yticks([])

plt.tight_layout()
#plt.savefig('../../../Desktop/seq_replay_PV1069.pdf')

#%% Get skewness value in test data
actual_wake_skew_vals = skew(H_wake,axis=1)
actual_REM_skew_vals = skew(H_REM,axis=1)

shuffled_wake_skew_vals = np.zeros((numShuffles,K))
shuffled_REM_skew_vals = np.zeros((numShuffles,K))

for shuffle_i in tqdm(range(numShuffles)):
    shuffled_W = np.zeros(W.shape)
    for neuron in range(numNeurons):
        shuffled_W[neuron,:,:] = np.roll(W[neuron,:,:],shift=np.random.randint(W.shape[2]),axis=1)
    
    shuffled_H_wake = extract_H(shuffled_W,wake_data)
    shuffled_H_REM = extract_H(shuffled_W,full_REM_data)
    shuffled_wake_skew_vals[shuffle_i,:] = skew(shuffled_H_wake,axis=1)
    shuffled_REM_skew_vals[shuffle_i,:] = skew(shuffled_H_REM,axis=1)

#[p_value] = test_significance(testData,W,n_surrogates=1000)

#%% Plot histogram of shuffles
plt.figure(figsize=(.75,.75))
plt.hist(shuffled_wake_skew_vals[:,0], histtype='stepfilled',alpha=.5,color='k',label='Shuffled')
plt.plot([actual_wake_skew_vals[0],actual_wake_skew_vals[0]],[0,200],color='C0',label='Actual')
plt.title('Wake, S1')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../../Desktop/PV1069_awake_S1_hist.pdf')

plt.figure(figsize=(.75,.75))
plt.hist(shuffled_wake_skew_vals[:,1], histtype='stepfilled',alpha=.5,color='k',label='Shuffled')
plt.plot([actual_wake_skew_vals[1],actual_wake_skew_vals[1]],[0,200],color='C6',label='Actual')
plt.title('Wake, S2')
plt.xlabel('Sequence distribution\nskewness')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../../Desktop/PV1069_awake_S2_hist.pdf')

#%% Plot histogram of shuffles
plt.figure(figsize=(.75,.75))
plt.hist(shuffled_REM_skew_vals[:,0], histtype='stepfilled',alpha=.5,color='k',label='Shuffled')
plt.plot([actual_REM_skew_vals[0],actual_REM_skew_vals[0]],[0,200],color='C0',label='Actual')
plt.title('REM, S1')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../../Desktop/PV1069_REM_S1_hist.pdf')

plt.figure(figsize=(.75,.75))
plt.hist(shuffled_REM_skew_vals[:,1], histtype='stepfilled',alpha=.5,color='k',label='Shuffled')
plt.plot([actual_REM_skew_vals[1],actual_REM_skew_vals[1]],[0,200],color='C6',label='Actual')
plt.title('REM, S2')
plt.xlabel('Sequence distribution\nskewness')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../../Desktop/PV1069_REM_S2_hist.pdf')

#%% Wake
print(f'S1 pvalue: p = {sum(shuffled_wake_skew_vals[:,0]>actual_wake_skew_vals[0])/numShuffles}')
print(f'S2 pvalue: p = {sum(shuffled_wake_skew_vals[:,1]>actual_wake_skew_vals[1])/numShuffles}')

#%% REM
print(f'S1 pvalue: p = {sum(shuffled_REM_skew_vals[:,0]>actual_REM_skew_vals[0])/numShuffles}')
print(f'S2 pvalue: p = {sum(shuffled_REM_skew_vals[:,1]>actual_REM_skew_vals[1])/numShuffles}')
# %%
