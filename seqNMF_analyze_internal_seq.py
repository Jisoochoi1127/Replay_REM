#%% Imports
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('plot_style.mplstyle')

import pingouin as pg
import h5py

#%%
results_dir = '../../output_REM/'
resultsList=os.listdir(results_dir)

# %% Q1: is there any structure activity during REM?
data_list = []
for file_name in resultsList:
    if file_name.startswith('seqResults_') and file_name.endswith('.h5'): # Only include seqResults, not replay results
        h5_file = h5py.File(os.path.join(results_dir,file_name))
        data_list.append( #This will create one list entry per cell
                {
                    'state':h5_file['state'][()].decode("utf-8"),
                    'condition':h5_file['condition'][()].decode("utf-8"),
                    'mouse':h5_file['mouse'][()].decode("utf-8"),
                    'S1_numSeqs':h5_file['S1_numSeqs'][()],
                    'S2_numSeqs':h5_file['S2_numSeqs'][()],
                    'S1_score':h5_file['S1_score'][()],
                    'S2_score':h5_file['S2_score'][()],
                    'S1_pvalue':h5_file['S1_pvalue'][()],
                    'S2_pvalue':h5_file['S2_pvalue'][()],
                }
            )

        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)

df_stats = df
df_stats.loc[df_stats['S1_pvalue']>0.05,'S1_numSeqs']=0
df_stats.loc[df_stats['S2_pvalue']>0.05,'S2_numSeqs']=0

#%%
df_numSeqs=df_stats.melt(id_vars=['state','mouse', 'condition'],value_name='numSeqs',value_vars=['S1_numSeqs', 'S2_numSeqs'],var_name='seqType')
df_seqScore=df.melt(id_vars=['state','mouse', 'condition'],value_name='seqScore',value_vars=['S1_score', 'S2_score'],var_name='seqType')

# %%
sns.barplot(
    data=df_numSeqs,
    y='numSeqs',
    x='condition',
    hue='state',
    order=['LTD1', 'LTD5', 'HATD1', 'HATD5'],
    hue_order=['REMpre','wake','REMpost'],
    palette=['C1','C4','C2'],
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_numSeqs,
    y='numSeqs',
    x='condition',
    hue='state',
    palette=['C1','C4','C2'],
    # color='k',
    order=['LTD1', 'LTD5', 'HATD1', 'HATD5'],
    hue_order=['REMpre','wake','REMpost'],
    size=1,
    dodge=True,
    legend=False
)
plt.xticks([0,1,2,3],['Novel', 'Fam.', 'Anxiety D1', 'Anxiety D5'],rotation=90)
plt.ylabel('Num. significant \nsequences')
plt.xlabel('')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/numSigSequences.pdf')

#%%STATS

#%% Sequence score
sns.barplot(
    data=df_seqScore,
    y='seqScore',
    x='condition',
    hue='state',
    order=['LTD1', 'LTD5', 'HATD1', 'HATD5'],
    hue_order=['REMpre','wake','REMpost'],
    palette=['C1','C4','C2'],
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_seqScore,
    y='seqScore',
    x='condition',
    hue='state',
    palette=['C1','C4','C2'],
    order=['LTD1', 'LTD5', 'HATD1', 'HATD5'],
    hue_order=['REMpre','wake','REMpost'],
    size=1,
    dodge=True,
    legend=False
)
plt.xticks([0,1,2,3],['Novel', 'Fam.', 'Anxiety D1', 'Anxiety D5'],rotation=90)
plt.ylabel('Num. significant \nsequences')
plt.xlabel('')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/seqScores.pdf')

#%% STATS



