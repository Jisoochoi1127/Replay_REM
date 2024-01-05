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
        data['binaryData'] = data['rawData']

    elif 'post' in state:
        filename='all_binary_post_REM'
        data = open_file(path, filename)
        data['binaryData'] = data['rawData']

    else: # Then, must be wake data in old mat format
        filename='ms'
        data = open_file(path, filename)
        data['binaryData'], _ = binarize_ca_traces(data['rawData'], 2, params['samplingFrequency'])
    
    return data

#%% TEMP TEST
mouse='pv1060'
condition='LTD5'
state='REMpre'
data = load_data(mouse, condition, state, params)

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
import matplotlib
#cmap = matplotlib.cm.get_cmap('nipy_spectral')
cmap = matplotlib.cm.get_cmap('viridis')
#numNeurons=data['SFPs'].shape[0]
numNeurons=100
threshold=.2

plt.figure(figsize=(3,1))
plt.subplot(131)
plt.imshow(REMpre_data['binaryData'].T, cmap='GnBu', aspect='auto', interpolation='none')
plt.xlim(0,1800)
plt.xticks([0,900,1800],[0,30,60])
plt.ylim(0,numNeurons)
#plt.xlabel('Time (s)')
plt.ylabel('Neuron #')
plt.title('REM pre')

plt.subplot(132)
plt.imshow(wake_data['binaryData'].T, cmap='GnBu', aspect='auto', interpolation='none')

plt.xlim(0,1800)
plt.xticks([0,900,1800],[0,30,60])
plt.ylim(0,numNeurons)
plt.yticks([])
plt.xlabel('Time (s)')
plt.ylabel('')
plt.title('Wakefulness')

plt.subplot(133)
plt.imshow(REMpost_data['binaryData'].T, cmap='GnBu', aspect='auto', interpolation='none')

plt.xlim(0,1800)
plt.xticks([0,900,1800],[0,30,60])
plt.ylim(0,numNeurons)
#plt.xlabel('Time (s)')
plt.yticks([])
plt.ylabel('')
plt.title('REM post')

plt.savefig('../../output_REM/calcium_transients.pdf')

#%% Plot anxiety track avoidance
HAT_data={}
path = '../../datasets/REM_data/pv1060/HATD5'

f = sio.loadmat(path + '/behav.mat')
HAT_data.update(
    { # Note that older files do not have background/tones
    'position':f['behav']['position'][0][0],
    'behavTime':f['behav']['time'][0][0].T[0]/1000,
    }
    )

HAT_data['LT_position'] = -1*(HAT_data['position'][:,0]-100)

plt.figure(figsize=(1.5,.5))
plt.subplot(1,2,1)
plt.plot(HAT_data['behavTime'],HAT_data['LT_position'],color='C4')
plt.xticks([0,300])
plt.xlabel('Time (s)')
plt.ylabel('Position\n(cm)')

plt.subplot(1,4,3)
sns.kdeplot(y=HAT_data['LT_position'],
            color='C4'
             )
plt.xticks([0,.025],[0,.025])
plt.yticks([])
plt.savefig('../../output_REM/anxiety_behavior.pdf')

#%% Plot anxiety track avoidance
normal_data={}
path = '../../datasets/REM_data/pv1060/LTD5'

f = sio.loadmat(path + '/behav.mat')
normal_data.update(
    { # Note that older files do not have background/tones
    'position':f['behav']['position'][0][0],
    'behavTime':f['behav']['time'][0][0].T[0]/1000,
    }
    )

normal_data['LT_position'] = -1*(normal_data['position'][:,0]-100)

plt.figure(figsize=(1.5,.5))
plt.subplot(1,2,1)
plt.plot(normal_data['behavTime'],normal_data['LT_position'])
plt.xticks([0,300])
plt.xlabel('Time (s)')
plt.ylabel('Position\n(cm)')

plt.subplot(1,4,3)
sns.kdeplot(y=normal_data['LT_position'],
             )
plt.xticks([0,.025],[0,.025])
plt.yticks([])
plt.savefig('../../output_REM/normal_behavior.pdf')
