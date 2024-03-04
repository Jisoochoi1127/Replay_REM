# %% Imports
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import yaml
import h5py
from scipy.ndimage import zoom
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

# %%
W_LT, H_LT, _, _, _ = seqnmf(
    data_LT['binaryData'][:,selected_neurons].T,
    K=params['K'],
    L=params['L'],
    Lambda=params['Lambda'],
    max_iter=10
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
mouse_list = ["pv1043", "pv1060", "pv1069", "pv1254"]

expansion_scores = np.zeros((len(expansionFactorList),len(mouse_list)))
expansion_numSeqs = np.zeros((len(expansionFactorList),len(mouse_list)))
params['maxIters'] = 10

for mouse_i, mouse in enumerate(mouse_list):
    data_LT = load_data(mouse=mouse, condition=condition, state="wake", params=params)
    data_REMpost = load_data(
        mouse=mouse, condition=condition, state="REMpost", params=params
    )

    with h5py.File(
        os.path.join(
            params["path_to_output"],
            "neuron_selection",
            f"selected_neurons_{condition}_{mouse}.h5",
        ),
        "r",
    ) as f:
        selected_neurons = f["place_cells"][()]

    for factor_i, factor in enumerate(tqdm(expansionFactorList)):
        params['expansionFactor'] = factor
        seqReplay_scores, seqReplay_pvalues, seqReplay_locs, _ = extract_seqReplay_score(data_LT['binaryData'][:,selected_neurons],
                                                                                    data_REMpost['binaryData'][:,selected_neurons],
                                                                                    params)

        expansion_scores[factor_i,mouse_i] = (seqReplay_scores[0]+seqReplay_scores[1])/2
        expansion_numSeqs[factor_i,mouse_i] = (len(seqReplay_locs[0])+len(seqReplay_locs[1]))/2

    with h5py.File(os.path.join(params["path_to_output"], f"{mouse}_sequence_expansion.h5"), "w") as f:
        f.create_dataset("expansion_scores", data=expansion_scores)
        f.create_dataset("expansion_numSeqs", data=expansion_numSeqs)
    
# %% Plot
plt.plot(expansionFactorList, expansion_scores)
plt.xlabel("Scaling factor")
plt.ylabel("Seq. score")
plt.savefig("../../output_REM/expansion_scores.pdf")

# %% Plot
plt.plot(expansionFactorList, expansion_numSeqs)
plt.xlabel("Scaling factor")
plt.ylabel("Num. sig.\nsequences")
plt.savefig("../../output_REM/expansion_numSeq.pdf")
# %%
