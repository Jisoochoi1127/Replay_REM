#%% Imports
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('plot_style.mplstyle')

import pingouin as pg

#%%
df = pd.read_csv('../../output_REM/seqScores.csv')
# %%
sns.barplot(
    data=df,
    y='max_seq_zscore',
    x='condition',
    hue='state',
    palette=['C0','C6','C1'],
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df,
    y='max_seq_zscore',
    x='condition',
    hue='state',
    palette=['C0','C6','C1'],
    size=2,
    dodge=True,
    legend=False
)
plt.xticks([0,1,2],['Novel', 'Fam.', 'Anxiety'])
plt.ylabel('Seq. score (z)')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/seq_zscore.png')

#%%STATS


#%%
sns.barplot(
    data=df,
    y='max_seq_score',
    x='condition',
    hue='state',
    palette=['C0','C6','C1'],
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df,
    y='max_seq_score',
    x='condition',
    hue='state',
    palette=['C0','C6','C1'],
    size=2,
    dodge=True,
    legend=False
)
plt.ylabel('Seq. score')
plt.xticks([0,1,2],['Novel', 'Fam.', 'Anxiety'])
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/seq_score.png')
#%% STATS


#%% Load replay data
df_replay = pd.read_csv('../../output_REM/seqReplayScores.csv')

# %%
sns.heatmap(data=df_replay.query("condition=='LTD1'").pivot_table(index='state_ref', 
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
plt.title('Novelty')
plt.xlabel('Reference')
plt.ylabel('Target')
plt.savefig('../../output_REM/replay_novelty.png')
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
