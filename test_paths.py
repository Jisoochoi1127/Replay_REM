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
from tqdm import tqdm
import itertools
from utils.helperFunctions import load_data, open_file

#%%
with open('params.yaml','r') as file:
    params = yaml.full_load(file)

#%% Define conditions, dataset
states_list = ['REMpre', 'wake', 'REMpost']
condition_list = ["LTD1", "LTD5", "HATD1", "HATD5", "HATDSwitch"]
mouse_list = ["pv1043", "pv1060", "pv1069", "pv1191", "pv1192", "pv1252", "pv1254"]

#%% TEMP TEST
for condition, mouse, state in tqdm(list(itertools.product(condition_list,
                                                           mouse_list,
                                                           states_list)),
                                                           total=len(condition_list)*len(mouse_list)*len(states_list)):
    if os.path.exists(os.path.join(params['path_to_dataset'],mouse,condition)):
        print(f'processing: {mouse},{condition},{state}')
        data = load_data(mouse, condition, state, params)
# %%
