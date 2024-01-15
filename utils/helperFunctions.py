import numpy as np
from seqnmf import seqnmf
import os
from scipy.stats import skew
import h5py
import scipy.io as sio
from scipy.signal import find_peaks
from pycaan.functions.signal_processing import preprocess_data, binarize_ca_traces
from pycaan.functions.dataloaders import load_data as pycaan_load

def open_file(path, filename):
    data={}
    try:
        f = h5py.File(os.path.join(path,filename)+'.mat', 'r')
        data.update({'data':f[filename][()].T})
    except:
        f = sio.loadmat(os.path.join(path,filename)+'.mat')
        try:
            data.update({'data':f[filename]['RawTraces'][0][0]})
        except:
            data.update({'data':f[filename]})
    return data

def load_data(mouse, condition, state, params):
    path=os.path.join(params['path_to_dataset'], mouse, condition)
    if 'pre' in state:
        filename='all_binary_pre_REM'
        data = open_file(path, filename)
        data['binaryData'] = data['data']

    elif 'post' in state:
        filename='all_binary_post_REM'
        data = open_file(path, filename)
        data['binaryData'] = data['data']

    elif state=='wake': # Then, must be wake data in old mat format
        data = preprocess_data(pycaan_load(path), params)
        data['binaryData'], _ = binarize_ca_traces(data['RawData'], 2, params['samplingFrequency'])

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