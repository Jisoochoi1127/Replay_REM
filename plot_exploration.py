#%%
import numpy as np
import os
from tqdm import tqdm
import itertools

import yaml
import h5py
from utils.helperFunctions import load_data
import matplotlib.pyplot as plt
plt.style.use('plot_style.mplstyle')

#%% Load parameters
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

np.random.seed(params['seed'])

#%% Define conditions, dataset
condition_list = ['LTD1','LTD5','HATD1','HATD5']
mouse_list = ['pv1043','pv1060', 'pv1069', 'pv1191', 'pv1192', 'pv1252', 'pv1254']
# %%
ct=1
plt.figure(figsize=(len(condition_list),len(mouse_list)))
for condition, mouse in tqdm(list(itertools.product(condition_list,
                                                           mouse_list
                                                           )),
                                                           total=len(condition_list)*len(mouse_list)):
    # try:
        
    plt.subplot(len(mouse_list),len(condition_list),ct)

    data = load_data(mouse=mouse, condition=condition, state='wake',params=params)
    plt.plot(data['position'])
        
    # except:
        # print('Missing session')
    
    ct+=1
# %%
