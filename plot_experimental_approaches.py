#%%
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('plot_style.mplstyle')

import pingouin as pg
import h5py
import yaml
import scipy.io as sio
from pycaan.functions.signal_processing import binarize_ca_traces

#%%
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
        'rawData':f[filename][()].T,
        }
        )
    except:
        f = sio.loadmat(path+'/ms.mat')
        data.update(
        {
        'rawData':f['ms']['RawTraces'][0][0],
        'SFPs':f['ms']['SFPs'][0][0],
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
    
    return data

#%% First, plot SFPs
mouse='pv1060'
condition='LTD1'
state='wake'
data = load_data(mouse, condition, state, params)

#%%
plt.figure()
plt.imshow(np.max(data['SFPs'],axis=2),
           vmin=0,
           vmax=10,
           cmap='YlGnBu_r')
plt.axis('off')
plt.savefig('../../output_REM/SFPs.pdf')
#%% Then plot transients during each state
mouse='pv1060'
condition='LTD1'
state='wake'
wake_data = load_data(mouse, condition, state, params)
state='REMpre'
REMpre_data = load_data(mouse, condition, state, params)
state='REMpost'
REMpost_data = load_data(mouse, condition, state, params)

#%%



# %%
