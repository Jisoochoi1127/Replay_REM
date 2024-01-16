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
from tqdm import tqdm


#%%
results_dir = '../../output_REM/'
resultsList=os.listdir(results_dir)

data_list = []

tuning_curves = []

for file_name in tqdm(resultsList):
    if file_name.startswith('tuning_') and file_name.endswith('.h5') and 'pv1254' not in file_name: # Only include seqResults, not replay results
        h5_file = h5py.File(os.path.join(results_dir,file_name))
        for i in range(len(h5_file['info'][()])):
            data_list.append( #This will create one list entry per cell
                {
                    'cell_ID': i,
                    'mouse': h5_file['mouse'][()].decode("utf-8"),
                    'condition': h5_file['condition'][()].decode("utf-8"),
                    'info': h5_file['info'][i],
                    'p_value': h5_file['p_value'][i],
                    'peak_loc': h5_file['peak_loc'][i][0],
                    'peak_val': h5_file['peak_val'][i],
                    'marginal_likelihood': h5_file['marginal_likelihood'][i],
                    'trans_prob': h5_file['trans_prob'][i],
                    'total_distance_travelled': h5_file['total_distance_travelled'][()]
                }
            )
            tuning_curves.append(h5_file['tuning_curves'][i])

        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)
tuning_curves = np.array(tuning_curves)

#%% plot place cells on linear track
sns.histplot(data=df.query("p_value<0.05"),
             x='peak_loc',
             y='condition',
             cbar=True,
             cbar_kws={'label':'Number of centroids'},
             cmap='magma')
plt.title('Place cell centroid locations')
plt.xlabel('Location (cm)')
# %%
sns.histplot(data=df.query("p_value<0.05 and condition=='HATD1' or condition=='HATD5'"),
             x='peak_loc',
             y='peak_val',
             cbar=True,
             cbar_kws={'label':'Number of centroids'})
# %%
sns.histplot(data=df.query("p_value<0.05"),
             x='info'
)

#%% 
plt.scatter(data=df,
            x='info',
            y='p_value')
plt.plot([0,.5],[.05,.05],'C4:')
plt.xlabel('Info.')
plt.ylabel('p value')


#%%