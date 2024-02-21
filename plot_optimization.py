# %% Imports
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import h5py
import yaml
from tqdm import tqdm

plt.style.use("plot_style.mplstyle")

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %%
mouseList = ['pv1060','pv1069','pv1254', 'pv1043']
L_list = [25, 50, 75, 125, 150, 175, 200, 225, 250, 300]

optimal_L_scores = []
for mouseName in mouseList:
    with h5py.File(f"../../output_REM/{mouseName}_optimal_L.h5") as f:
        for i_L, L in enumerate(L_list):
            optimal_L_scores.append(
                {
                    "mouse": mouseName,
                    "L": L,
                    "S1_numSeq": f['L_scores'][i_L][0],
                    "S2_numSeq": f['L_scores'][i_L][1]
                }
            )

df = pd.DataFrame(optimal_L_scores)

df_melted = df.melt(
    id_vars=['mouse', 'L'],
    value_vars=['S1_numSeq', 'S2_numSeq'],
    var_name='SeqType',
    value_name='numSeq'
)
#%% Plot results
sns.lineplot(
    data=df_melted,
    x='L',
    y='numSeq',
    errorbar='se',
    color='C1'
)
plt.title('Optimal L\n\
          n = 4 mice, 2 seq.')
plt.xlabel('$L$\n(frame number)')
plt.ylabel("Num. sig.\nsequences")
plt.savefig("../../output_REM/optimization_L.pdf")

# %% Optimize lambda
mouseList = ['pv1060','pv1069','pv1254', 'pv1043']
lambda_list = [10e-6,10e-5,10e-4,10e-3,10e-2,10e-1,1,10,10e2,10e3,10e4]

optimal_lambda_scores = []
for mouseName in mouseList:
    with h5py.File(f"../../output_REM/{mouseName}_optimal_lambda.h5") as f:
        for i_L, L in enumerate(L_list):
            optimal_lambda_scores.append(
                {
                    "mouse": mouseName,
                    "lambda": L,
                    "S1_numSeq": f['L_scores'][i_L][0], #TODO rename to lambda
                    "S2_numSeq": f['L_scores'][i_L][1] #TODO rename to lambda
                }
            )

df = pd.DataFrame(optimal_lambda_scores)

df_melted = df.melt(
    id_vars=['mouse', 'lambda'],
    value_vars=['S1_numSeq', 'S2_numSeq'],
    var_name='SeqType',
    value_name='numSeq'
)

#%% Plot results
sns.lineplot(
    data=df_melted,
    x='lambda',
    y='numSeq',
    errorbar='se',
    color='C1'
)
plt.title('Optimal $\lambda$\n\
          n = 4 mice, 2 seq.')
plt.xlabel('$lambda$')
plt.ylabel("Num. sig.\nsequences")
plt.savefig("../../output_REM/optimization_lambda.pdf")
# %%
