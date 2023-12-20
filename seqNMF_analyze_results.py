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
df_stats=df_stats.melt(id_vars=['state','mouse', 'condition'],value_name='numSeqs',value_vars=['S1_numSeqs', 'S2_numSeqs'],var_name='seqType')
# %%
sns.barplot(
    data=df_stats,
    y='numSeqs',
    x='condition',
    hue='state',
    hue_order=['REMpre','wake','REMpost'],
    palette=['C1','C4','C2'],
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_stats,
    y='numSeqs',
    x='condition',
    hue='state',
    palette=['C1','C4','C2'],
    hue_order=['REMpre','wake','REMpost'],
    size=2,
    dodge=True,
    legend=False
)
plt.xticks([0,1,2],['Novel', 'Fam.', 'Anxiety'])
plt.ylabel('Num. significant \nsequences')
plt.xlabel('')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/numSigSequences.pdf')

#%%STATS


#%% Replay analysis
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


df_replay_stats = df_replay
df_replay_stats.loc[df_stats['S1_pvalue']>0.05,'S1_numSeqs']=0
df_replay_stats.loc[df_stats['S2_pvalue']>0.05,'S2_numSeqs']=0

#%%
df_replay_stats=df_replay_stats.melt(id_vars=['state_ref','state_pred','mouse', 'condition'],value_name='numSeqs',value_vars=['S1_numSeqs', 'S2_numSeqs'],var_name='seqType')
#%%
sns.heatmap(data=df_replay_stats.query("condition=='LTD1' and seqType=='S2_numSeqs'").pivot_table(index='state_ref', 
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
plt.savefig('../../output_REM/num_replayed_seqs_novelty.png')
# %%
sns.heatmap(data=df_replay.query("condition=='LTD5'").pivot_table(index='state_ref', 
                                columns='state_pred', 
                                values='max_seqReplay_zscore'),
                                cmap='magma',
                                vmin=2,
                                vmax=4,
                                rasterized=True,
                                cbar_kws={'label':'Seq. replay\nscore (z)'}
)
#plt.xticks([])
#plt.yticks([])
plt.title('Familiarity')
plt.xlabel('Reference')
plt.ylabel('Target')
plt.savefig('../../output_REM/replay_familiar.png')
# %%
sns.heatmap(data=df_replay.query("condition=='HATD5'").pivot_table(index='state_ref', 
                                columns='state_pred', 
                                values='max_seqReplay_zscore'),
                                cmap='magma',
                                vmin=2,
                                vmax=4,
                                rasterized=True,
                                cbar_kws={'label':'Seq. replay\nscore (z)'}
)
#plt.xticks([])
#plt.yticks([])
plt.title('Anxiety')
plt.xlabel('Reference')
plt.ylabel('Target')
plt.savefig('../../output_REM/replay_anxiety.png')
# %% Pre-play only
sns.barplot(
    data=df_replay.query("state_ref=='REMpre' and state_pred != 'REMpre'"),
    y='max_seqReplay_zscore',
    x='condition',
    hue='state_pred',
    palette=['C0','C6','C1'],
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_replay.query("state_ref=='REMpre' and state_pred != 'REMpre'"),
    y='max_seqReplay_zscore',
    x='condition',
    hue='state_pred',
    palette=['C0','C6','C1'],
    size=2,
    dodge=True,
    legend=False
)
plt.plot([-0.5,2.5],[2,2],'k:')
plt.xticks([0,1,2],['Novel', 'Fam.', 'Anxiety'])
plt.ylabel('Seq. preplay\nscore (z)')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/preplay_by_condition.png')

# %%
sns.barplot(
    data=df_replay.query("state_ref=='wake' and state_pred != 'wake'"),
    y='max_seqReplay_zscore',
    x='condition',
    hue='state_pred',
    palette=['C0','C6','C1'],
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_replay.query("state_ref=='wake' and state_pred != 'wake'"),
    y='max_seqReplay_zscore',
    x='condition',
    hue='state_pred',
    palette=['C0','C6','C1'],
    size=2,
    dodge=True,
    legend=False
)
plt.plot([-0.5,2.5],[2,2],'k:')
plt.xticks([0,1,2],['Novel', 'Fam.', 'Anxiety'])
plt.ylabel('Seq. replay\nscore (z)')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/replay_by_condition.png')

# %%
