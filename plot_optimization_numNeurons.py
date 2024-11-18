# %% Imports
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import h5py
import yaml

plt.style.use("plot_style.mplstyle")

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

# %%
condition = "LTD1"
mouse_ID = "pv1069" # "pv1254", "pv1043", "pv1060"
numNeurons_list = [8, 16, 32, 64, 128, 256, 512]

data_dict = []
with h5py.File(os.path.join(params["path_to_output"], f"{mouse_ID}_optimal_numNeurons.h5"), "r") as f:
    for i, numNeurons in enumerate(numNeurons_list):
        data_dict.append(
        {
            "num_neurons": numNeurons,
            "num_S1":f["numNeurons_scores"][i,0],
            "num_S2":f["numNeurons_scores"][i,1],
        }
    )
    
df = pd.DataFrame(data_dict)

#%%
df_melted = df.melt(
    id_vars=['num_neurons'],
    value_vars=['num_S1', 'num_S2'],
    var_name='SeqType',
    value_name='numSeq'
)
#%% Plot results
sns.lineplot(
    data=df_melted,
    x='num_neurons',
    y='numSeq',
    errorbar='se',
    color='C1'
)
plt.plot([params['numNeurons'],params['numNeurons']],[0,48],'k:')
plt.xlabel('Num. neurons')
plt.ylabel("Num. sig.\nsequences")
plt.savefig("../../output_REM/optimization_numNeurons.pdf")

# %% Optimize lambda
mouseList = ['pv1060','pv1069','pv1254', 'pv1043']
lambda_list = [10e-6,10e-5,10e-4,10e-3,10e-2,10e-1,1,10,10e2,10e3,10e4]

optimal_lambda_scores = []
for mouseName in mouseList:
    with h5py.File(f"../../output_REM/{mouseName}_optimal_lambda.h5") as f:
        for i_L, Lambda in enumerate(lambda_list):
            optimal_lambda_scores.append(
                {
                    "mouse": mouseName,
                    "lambda": Lambda,
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
plt.xscale('log')
plt.title('Optimal $\lambda$\n\
          n = 4 mice, 2 seq.')
plt.xlabel('$lambda$')
plt.ylabel("Num. sig.\nsequences")
plt.savefig("../../output_REM/optimization_lambda.pdf")
# %%
