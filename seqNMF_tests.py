#%% Imports
import numpy as np
from seqnmf import seqnmf
#from utils.utils import test_significance
import os
from scipy.signal import find_peaks
from scipy.stats import skew
from tqdm import tqdm

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

#%% Load awake data
data={}
path = '../../datasets/calcium_imaging/CA1/M246/M246_LT_6'
f = h5py.File(path + '/ms.mat','r')
data.update({
            'caTime':np.array(f.get('ms/time'))[0]/1000, # convert ms->s
            'rawData':np.array(f.get('ms/RawTraces')).T})
data['task'] = 'LT'
params['z_threshold'] = 2
params['sampling_frequency'] = 30

LT_data, _ = binarize_ca_traces(data['rawData'], 2, 30)
LT_data = LT_data[:,0:params['numNeurons']]

#%% Establishing train/test portions
numFrames, numNeurons = LT_data.shape
# For awake data, further refine into running epochs
# Define train/test datasets
trainingFrames = np.zeros(len(LT_data), dtype=bool)
trainingFrames[0:int(params['train_test_ratio']*len(LT_data))] = True
testingFrames = ~trainingFrames

# Split dataset in two
train_data = LT_data[trainingFrames]
test_data = LT_data[testingFrames]

#%% First, extract putative sequences during REMpre
W_train, H_train, _, _, _ = seqnmf(
    train_data.T,
    K=params['K'],
    L=params['L'],
    Lambda=params['Lambda'],
    max_iter=params['maxIters']
    )

#%% Find what each neuron belongs to what sequence
sorting_index = np.argsort(np.argmax(np.concatenate((W_train[:,0,:],W_train[:,1,:]),axis=1),axis=1))

#%% Custom function to extract H given X and W
def extract_H(W,data):
    H = np.zeros((W.shape[1], data.shape[0]),dtype='float')
    for l in range(W.shape[2]):
        X_shifted = np.roll(data,shift=1-l,axis=0)
        H = H + (X_shifted @ W[:, :, l]).T #Matrix multiplication
    return H

#%% extract H on test data
H_test = extract_H(W_train, test_data)

#%% Plot resulting sequences
plt.imshow(test_data[:,sorting_index].T,aspect='auto',cmap='magma',vmax=.5)

#%%
# Shuffle
shuffled_test_H = np.zeros((params['numShuffles'], params['K'],len(test_data)))
for shuffle_i in tqdm(range(params['numShuffles'])):
    shuffled_W = np.zeros(W_train.shape)
    for neuron in range(params['numNeurons']):
        shuffled_W[neuron,:,:] = np.roll(W_train[neuron,:,:],shift=np.random.randint(W_train.shape[2]),axis=1)
    
    shuffled_test_H[shuffle_i,:,:] = extract_H(shuffled_W,test_data)

#%%
#%% Compute p-value
# Define confidence as 1-p_value
H_test_pvalue = np.zeros((params['K'],len(test_data)))
zscored_H = np.zeros((params['K'],len(test_data)))
H_test_confidence = np.zeros((params['K'],len(test_data)))

for k in range(params['K']):
    zscored_H[k,:] = (H_test[k,:]-np.mean(H_test[k,:]))/np.std(H_test[k,:])
    for i in range(len(test_data)):
        H_test_pvalue[k,i] = sum(shuffled_test_H[:,k,i]>H_test[k,i])/params['numShuffles']
    
    H_test_confidence[k,:] = 1-H_test_pvalue[k,:]
    H_test_confidence[k,zscored_H[k,:]<2] = 0

#%% Find peaks


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
