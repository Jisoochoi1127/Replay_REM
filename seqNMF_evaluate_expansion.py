# %% Imports
import numpy as np
import os
import matplotlib.pyplot as plt
import yaml
import h5py
from scipy.ndimage.interpolation import zoom
from seqnmf import seqnmf
from utils.helperFunctions import load_data, extract_seqReplay_score
from tqdm import tqdm

plt.style.use("plot_style.mplstyle")

# %% Params
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %% Load LT and REM data
mouse = 'pv1060'
condition = 'LTD1'
data_LT = load_data(mouse=mouse,
                    condition=condition,
                    state='wake',
                    params=params)
data_REMpost = load_data(mouse=mouse,
                    condition=condition,
                    state='REMpost',
                    params=params)

with h5py.File(os.path.join(params['path_to_output'],"neuron_selection", f'selected_neurons_{condition}_{mouse}.h5'),'r') as f:
    selected_neurons = f['place_cells'][()]

data_LT['binaryData'][:,selected_neurons]

# %% Create ~20 different compressed/streched templates by interpolation
params['maxIter'] = 10 # Temp
W_LT, H_LT, _, _, _ = seqnmf(
    data_LT['binaryData'][:,selected_neurons].T,
    K=params['K'],
    L=params['L'],
    Lambda=params['Lambda'],
    max_iter=params['maxIter']
    )

#%% Test temp
sortingIndex = np.argsort(np.argmax(W_LT[:,0,:],axis=1),axis=0)

shrunk = zoom(W_LT, (1,1,0.5))
expanded = zoom(W_LT, (1,1,2))
plt.figure(figsize=(2,1))
plt.subplot(182)
plt.title('Half')
plt.imshow(shrunk[sortingIndex,0,:], interpolation='none', aspect='auto', vmin=0,vmax=.5)
plt.axis('off')

plt.subplot(142)
plt.title('Original')
plt.imshow(W_LT[sortingIndex,0,:],interpolation='none',aspect='auto',vmin=0,vmax=.5)
plt.axis('off')

plt.subplot(122)
plt.title('Double')
plt.imshow(expanded[sortingIndex,0,:],interpolation='none', aspect='auto', vmin=0,vmax=.5)
plt.axis('off')
plt.savefig("../../output_REM/seq_expansion.pdf")

# %% For each template, compute matrix multiplication to extract H
expansionFactorList = np.arange(.5,2.1,.1)
expansion_scores = np.zeros(len(expansionFactorList))
expansion_numSeqs = np.zeros(len(expansionFactorList))
params['maxIter'] = 10
for i, factor in enumerate(tqdm(expansionFactorList)):
    params['expansionFactopr'] = factor
    seqReplay_scores, seqReplay_pvalues, seqReplay_locs = extract_seqReplay_score(data_LT['binaryData'][:,selected_neurons], data_REMpost['binaryData'][:,selected_neurons], params)

    expansion_scores[i] = (seqReplay_scores[0]+seqReplay_scores[1])/2
    expansion_numSeqs[i] = (len(seqReplay_locs[0])+len(seqReplay_locs[1]))/2

# %% Compute