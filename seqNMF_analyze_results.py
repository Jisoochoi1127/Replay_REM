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
    errorbar='se'
)
plt.ylabel('Seq. score (z)')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)

#%%
sns.barplot(
    data=df,
    y='max_seq_score',
    x='condition',
    hue='state',
    palette=['C0','C6','C1'],
    errorbar='se'
)
plt.ylabel('Seq. score')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)

#%% STATS
