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

data_list = []
for file_name in resultsList:
    if file_name.startswith('seqReplayResults_') and file_name.endswith('.h5'): # Only include seqResults, not replay results
        h5_file = h5py.File(os.path.join(results_dir,file_name))
        data_list.append( #This will create one list entry per cell
                {
                    'state_ref':h5_file['state_ref'][()].decode("utf-8"),
                    'state_pred':h5_file['state_pred'][()].decode("utf-8"),
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

df_replay = pd.DataFrame(data_list)

#%%
df_replay_stats = df_replay
df_replay_stats.loc[df_replay_stats['S1_pvalue']>0.05,'S1_numSeqs']=0
df_replay_stats.loc[df_replay_stats['S2_pvalue']>0.05,'S2_numSeqs']=0
df_replay_stats=df_replay_stats.melt(id_vars=['state_ref','state_pred','mouse', 'condition'],value_name='numSeqs',value_vars=['S1_numSeqs', 'S2_numSeqs'],var_name='seqType')


#%% Preplay?
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_replay_stats.query("state_ref == 'REMpre' and \
                                state_pred == 'wake' and \
                               condition == 'LTD5'"),
    x='mouse',
    y='numSeqs',
    color='C1',
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_replay_stats.query("state_ref == 'REMpre' and \
                                state_pred == 'wake' and \
                               condition == 'LTD5'"),
    x='mouse',
    y='numSeqs',
    hue='seqType',
    # color='C0',
    palette=(['C3','C4']),
    size=2,
    #dodge=True,
    legend=True
)
#plt.xticks([0,1],['S1', 'S2'])
plt.ylabel('Num. significant \nsequences')
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/numSeqPreplay_LT.pdf')

#%% Effect of conditions
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_replay_stats.query("state_ref == 'wake' and \
                                state_pred == 'REMpost' and \
                               condition != 'HATD1'"),
    x='condition',
    y='numSeqs',
    hue='seqType',
    palette=(['C0','C4']),
    color='C1',
    errorbar='se',
    capsize=.2,
)
sns.stripplot(
    data=df_replay_stats.query("state_ref == 'wake' and \
                                state_pred == 'REMpost' and \
                               condition != 'HATD1'"),
    x='condition',
    y='numSeqs',
    hue='seqType',
    # color='C0',
    palette=(['C0','C4']),
    size=2,
    dodge=True,
    legend=False
)
#plt.xticks([0,1],['S1', 'S2'])
plt.ylabel('Num. significant \nsequences')
plt.xlabel('')
plt.xticks([0,1,2],['Anxiety','Novelty','Familiarity'],rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/numSeqPreplay_conditions.pdf')


#%%
sns.heatmap(data=df_replay_stats.query("condition=='LTD1'").pivot_table(index='state_ref', 
                                columns='state_pred', 
                                values='numSeqs'),
                                cmap='magma',
                                rasterized=True,
                                cbar_kws={'label':'Num. replayed\nsequences'}
)
#plt.xticks([])
#plt.yticks([])
plt.title('Novelty')
plt.xlabel('Reference')
plt.ylabel('Target')
plt.savefig('../../output_REM/numSeqs_novelty_heatmap.pdf')
# %%
sns.heatmap(data=df_replay_stats.query("condition=='LTD5'").pivot_table(index='state_ref', 
                                columns='state_pred', 
                                values='numSeqs'),
                                cmap='magma',
                                rasterized=True,
                                cbar_kws={'label':'Num. replayed\nsequences'}
)
#plt.xticks([])
#plt.yticks([])
plt.title('Familiarity')
plt.xlabel('Reference')
plt.ylabel('Target')
plt.savefig('../../output_REM/numSeqs_familiar_heatmap.pdf')
# %%
sns.heatmap(data=df_replay_stats.query("condition=='HATD1'").pivot_table(index='state_ref', 
                                columns='state_pred', 
                                values='numSeqs'),
                                cmap='magma',
                                rasterized=True,
                                cbar_kws={'label':'Num. replayed\nsequences'}
)
#plt.xticks([])
#plt.yticks([])
plt.title('Novel anxiety')
plt.xlabel('Reference')
plt.ylabel('Target')
plt.savefig('../../output_REM/numSeqs_novelAnxiety_heatmap.pdf')

# %%
sns.heatmap(data=df_replay_stats.query("condition=='HATD5'").pivot_table(index='state_ref', 
                                columns='state_pred', 
                                values='numSeqs'),
                                cmap='magma',
                                rasterized=True,
                                cbar_kws={'label':'Num. replayed\nsequences'}
)
#plt.xticks([])
#plt.yticks([])
plt.title('Familiar anxiety')
plt.xlabel('Reference')
plt.ylabel('Target')
plt.savefig('../../output_REM/numSeqs_familiarAnxiety_heatmap.pdf')




# %% Pre-play only
plt.figure(figsize=(1,1))
sns.barplot(
    data=df_replay_stats.query("state_ref=='wake' and state_pred == 'REMpost'"),
    y='numSeqs',
    hue='condition',
    x='seqType',
    palette=['C4','C3','C1','C0'],
    errorbar='se',
    capsize=.2,
)
sns.stripplot(
    data=df_replay_stats.query("state_ref=='wake' and state_pred == 'REMpost'"),
    y='numSeqs',
    hue='condition',
    x='seqType',
    palette=['C4','C3','C1','C0'],
    size=1,
    dodge=True,
    legend=False
)

plt.xticks([0,1],['Seq. 1', 'Seq. 2'])
plt.xlabel('Sequence type')
plt.ylabel('Num. seq.\nreplayed')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/replay_by_condition_and_seqType_bars.pdf')