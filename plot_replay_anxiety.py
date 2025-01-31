# %% Imports
import os
import numpy as np
import pandas as pd
import seaborn as sns
import pingouin as pg
import matplotlib.pyplot as plt
import yaml
import h5py
from tqdm import tqdm
from utils.helperFunctions import load_data

plt.style.use("plot_style.mplstyle")

#%% Load params
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

#%% Assemblies
results_dir = params['path_to_output']+"/assembly"
resultsList = os.listdir(results_dir)

#%%
data_list = []
data_list_num_events = []

for file_name in tqdm(resultsList):
    if (
        file_name.startswith("assembly_")
        and file_name.endswith(".h5")
        and "pv1254" not in file_name # Exclude pv1254
    ):
        h5_file = h5py.File(os.path.join(results_dir, file_name))
        data = load_data(h5_file["mouse"][0].decode("utf-8"),
                         h5_file["condition"][0].decode("utf-8"),
                         'REMpost',
                         params
                         )
        numFrames = data['binaryData'].shape[0]
        recordingLength = numFrames/params['sampling_frequency']

        data_list_num_events.append(
            {
                "mouse": h5_file["mouse"][0].decode("utf-8"),
                "condition": h5_file["condition"][0].decode("utf-8"),
                "numSigEvents": np.sum(h5_file['post_rem_A_sig'][()]==1),
                "numAssemblies": np.max(h5_file['post_rem_A_ID']),
                "meanReactPerAssembly": np.mean(np.unique(h5_file['post_rem_A_ID'][()],return_counts=True)),
                "meanFreqPerAssembly": np.mean(np.unique(h5_file['post_rem_A_ID'][()],return_counts=True))/recordingLength
            }
        )

        for i in range(len(h5_file["post_rem_react_idx"][()])):
            data_list.append(  # This will create one list entry per cell
                {
                    "eventID": i,
                    "mouse": h5_file["mouse"][0].decode("utf-8"),
                    "condition": h5_file["condition"][0].decode("utf-8"),
                    "replayEventTime": h5_file['post_rem_react_idx'][i]/params['sampling_frequency'],
                    "replayEventID": h5_file['post_rem_A_ID'][i],
                    "replayEventStrength": h5_file['post_rem_A_str'][i],
                    "replayEventSignificance": h5_file['post_rem_A_sig'][i],
                    "peak_loc": h5_file['map_loc'][h5_file['post_rem_A_ID'][i]-1],
                    "wake_rate": h5_file['wake_rate'][h5_file['post_rem_A_ID'][i]-1],
                }
            )

        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)
df_numEvents = pd.DataFrame(data_list_num_events)

#%% Plot num assemblies
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_numEvents.query("condition=='LTD5' or condition=='HATD1'"),
    x='condition',
    y='numAssemblies',
    order=['LTD5', 'HATD1'],
    palette=(['C0','C4']),
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_numEvents.query("condition=='LTD5' or condition=='HATD1'"),
    color='gray',
    order=['LTD5', 'HATD1'],
    x='condition',
    y='numAssemblies',
    size=2
)
plt.xticks([0,1],['Control','Anxiety'],rotation=90)
plt.xlabel('')
plt.ylabel('Num. assemblies')
plt.savefig("../../output_REM/anxietyReplay_numAssemblies.pdf")

#%% Plot assembly frequency
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_numEvents.query("condition=='LTD5' or condition=='HATD1'"),
    x='condition',
    y='meanFreqPerAssembly',
    order=['LTD5', 'HATD1'],
    palette=(['C0','C4']),
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_numEvents.query("condition=='LTD5' or condition=='HATD1'"),
    color='gray',
    order=['LTD5', 'HATD1'],
    x='condition',
    y='meanFreqPerAssembly',
    size=2
)

plt.xticks([0,1],['Control','Anxiety'],rotation=90)
plt.xlabel('')
plt.ylabel('Mean reactivation\nfrequency (Hz)')
plt.savefig("../../output_REM/anxietyReplay_meanAssemblyFreq.pdf")

#%% Plot assembly prefered locs.
# plt.figure(figsize=(.75,1))
sns.histplot(
    data=df.query("condition=='LTD5' or condition=='HATD1'"),
    hue='condition',
    x='peak_loc',
    hue_order=['LTD5', 'HATD1'],
    palette=['C0','C4'],
    stat='density',
    legend=False,
    element='step'
)
plt.title('Distribution of reactivations')
plt.xlabel('Location on track (cm)')
plt.savefig("../../output_REM/anxietyReplay_assemblyPeakLoc.pdf")

#%% Separate rewards, open , closed locations

#%% Descriptives

#%% Stats


#%%




#%% Bayesian replay
if params['equalize_sampling']:
    results_dir = params['path_to_output']+"/equal_bayesian_replay"
else:
    results_dir = params['path_to_output']+"/bayesian_replay"
resultsList = os.listdir(results_dir)

data_list = []
num_event_list = []

for file_name in tqdm(resultsList):
    if (
        file_name.startswith("bayesian_replay_")
        and file_name.endswith(".h5")
        and "pv1254" not in file_name # Exclude pv1254
    ):
        h5_file = h5py.File(os.path.join(results_dir, file_name))
        num_event_list.append(
            {
                "mouse": h5_file["mouse"][()].decode("utf-8"),
                "condition": h5_file["condition"][()].decode("utf-8"),
                "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                "numReplayEvents": len(h5_file["replay_locs"][()])
            }
        )
        for i in range(len(h5_file["replay_locs"][()])):
            data_list.append(  # This will create one list entry per cell
                {
                    "eventID": i,
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "replayEventTime": h5_file['replay_locs'][i]/params['sampling_frequency'],
                    "replayEventScore": h5_file['replay_score'][i],
                    "replayEventJumpiness": h5_file['replay_jumpiness'][i],
                    "replayEventLength": h5_file['replay_length'][i],
                    "replayEventSlope": abs(h5_file['replay_slope'][i])
                }
            )
        
        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)
df_numEvents = pd.DataFrame(num_event_list)

#%% Plot number of events per condition
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_numEvents.query("Type=='replay' and condition=='HATD1' or condition=='LTD5'"),
    x='condition',
    y='numReplayEvents',
    order=['LTD5', 'HATD1'],
    palette=(['gray','C4']),
    #showfliers=False
    errorbar='se',
    capsize=.2
)

plt.ylabel('Num. replay\nevents')
plt.xticks([0,1],['Familiar','Anxiety'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/anxiety_bayesian_numEvents.pdf')

#%% DESCRIPTIVES (TODO: change var names to safe/anxiety)
mean_novelty = df_numEvents.query("condition=='LTD5' and Type=='replay'")['numReplayEvents'].mean()
SEM_novelty = df_numEvents.query("condition=='LTD5' and Type=='replay'")['numReplayEvents'].sem()
print(f'Novel: {mean_novelty} +/- {SEM_novelty}')

mean_familiar = df_numEvents.query("condition=='HATD1' and Type=='replay'")['numReplayEvents'].mean()
SEM_familiar = df_numEvents.query("condition=='HATD1' and Type=='replay'")['numReplayEvents'].sem()
print(f'Familiar: {mean_familiar} +/- {SEM_familiar}')

#%% STATS
pg.anova(
    data=df_numEvents.query("condition=='LTD5' or condition=='HATD1'"),
    dv='numReplayEvents',
    between='condition',

)

#%% Plot replay score as a function of novelty
plt.figure(figsize=(.75,1))
sns.boxenplot(
    data=df.query("Type=='replay' and condition=='LTD5' or condition=='HATD1'"),
    x='condition',
    y='replayEventScore',
    order=['LTD5', 'HATD1'],
    palette=(['gray','C4']),
    showfliers=False
    #errorbar='se',
    #capsize=.2
)

plt.ylabel('Replay score ($R^{2}$)')
plt.xticks([0,1],['Familiar','Anxiety'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/anxiety_bayesian_replay_scores.pdf')

#%% DESCRIPTIVES
mean_novelty = df.query("Type=='replay' and condition=='LTD5'")['replayEventScore'].mean()
SEM_novelty = df.query("Type=='replay' and condition=='LTD5'")['replayEventScore'].sem()
print(f'Novel: {mean_novelty} +/- {SEM_novelty}')

mean_familiar = df.query("Type=='replay' and condition=='HATD1'")['replayEventScore'].mean()
SEM_familiar = df.query("Type=='replay' and condition=='HATD1'")['replayEventScore'].sem()
print(f'Familiar: {mean_familiar} +/- {SEM_familiar}')

#%% STATS
pg.anova(
    data=df.query("Type=='replay' and condition=='LTD5' or condition=='HATD1'"),
    dv='replayEventScore',
    between='condition',
)

#%% Example posteriors
mouse = 'pv1069'
with h5py.File(os.path.join(params['path_to_output'],'equal_posterior_probs', 'posterior_probs_LTD5_pv1060.h5'), 'r') as f:
    familiar_posteriors = f[f'REMpost_posterior_probs'][()]

with h5py.File(os.path.join(params['path_to_output'],'equal_posterior_probs', 'posterior_probs_HATD1_pv1060.h5'), 'r') as f:
    anxiety_posteriors = f[f'REMpost_posterior_probs'][()]

plt.figure(figsize=(3.5,.75))
plt.subplot(141)
plt.imshow(
    familiar_posteriors.T,
    aspect='auto',
    interpolation='none',
    vmin=0.01,
    vmax=.065,
    cmap='Blues',
    rasterized=True
    )
plt.yticks([0,19,39],[100,50,0])
plt.ylabel('Location (cm)')
plt.xlim(0,4000)
plt.xticks([0,1800,3600],[0,1,2])
plt.xlabel('Time (min)')

plt.subplot(142)
plt.imshow(
    anxiety_posteriors.T,
    aspect='auto',
    interpolation='none',
    vmin=0.01,
    vmax=.065,
    cmap='Blues',
    rasterized=True
    )
plt.xlim(0,4000)
plt.yticks([0,19,39],['','',''])
plt.xticks([0,1800,3600],[0,1,2])

plt.subplot(185)
plt.imshow(np.zeros((1,1))*np.nan,
    vmin=0.01,
    vmax=.065,
    cmap='Blues',
           )
plt.axis('off')
plt.colorbar(label='Posterior prob.')

plt.subplot(188)
plt.plot(
    np.median(familiar_posteriors,axis=0),
    np.arange(40),
    color='C0'
    )
plt.plot(
    np.median(anxiety_posteriors,axis=0),
    np.arange(40),
    color='C4'
    )
plt.xlabel('Median posterior\nprobability')

plt.xlim(.022,.028)
plt.yticks([])
plt.xticks(rotation=90)
plt.savefig('../../output_REM/anxietyReplay_examplePosteriors.pdf')

# %% List posteriors results
if params['equalize_sampling']:
    results_dir = params['path_to_output']+'/equal_posterior_probs'
else:
    results_dir = params['path_to_output']+'/posterior_probs'
resultsList=os.listdir(results_dir)

# %% Load data
data_list = []
for file_name in tqdm(resultsList):
    if (file_name.startswith('posterior')
    and '1254' not in file_name
        ):
        with h5py.File(os.path.join(results_dir, file_name), 'r') as f:
            mouse = f['mouse'][()].decode("utf-8")
            condition = f['condition'][()].decode("utf-8")

            for state in ['REMpre', 'REMpost']:
                posterior_probs = f[f'{state}_posterior_probs'][()]
                avg_posteriors = np.mean(posterior_probs,axis=0)
                median_posteriors = np.median(posterior_probs,axis=0)
                max_posterior = np.argmax(posterior_probs,axis=1)

                for i, val in enumerate(avg_posteriors):
                    data_list.append(
                    {
                        'Mouse': mouse,
                        'Condition': condition,
                        'State': state,
                        'Location': (i+params['spatialBinSize']/2)*params['spatialBinSize'], # Convert to cm
                        'Average posterior': val,
                        'Median posterior': median_posteriors[i],
                        'MAP': max_posterior[i]
                    }
        )

df = pd.DataFrame(data_list)

# %% Median posterior
sns.lineplot(
    data=df.query("Condition =='LTD5' or Condition=='HATD1'"),
    x='Location',
    y='Median posterior',
    hue='Condition',
    hue_order=['LTD5', 'HATD1'],
    palette = ['C0','C4'],
    errorbar='se'
    )

plt.ylim(.023,.026)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/anxietyReplay_medianPosteriors.pdf')

#%% Bin open vs closed
df['Zone']=''
df.loc[df['Location']<=50,'Zone']='Closed'
df.loc[df['Location']>50,'Zone']='Open'

# %% Median posterior
sns.barplot(
    data=df.query("Condition =='LTD5' or Condition=='HATD1'"),
    x='Zone',
    y='Median posterior',
    hue='Condition',
    hue_order=['LTD5', 'HATD1'],
    palette = ['C0','C4'],
    errorbar='se',
    capsize=.2
    )

plt.ylim(.023,.026)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/anxietyReplay_medianPosteriors_zone_bars.pdf')

#%% STATS
pg.anova(
    data=df.query("Condition =='LTD5' or Condition=='HATD1'"),
    dv='Median posterior',
    between=['Zone','Condition'],
)

#%% POSTHOC
pg.pairwise_ttests(
    data=df.query("Condition =='LTD5' or Condition=='HATD1'"),
    dv='Median posterior',
    between=['Zone','Condition'],
)




# %% MAP
sns.lineplot(
    data=df.query("Condition =='LTD5' or Condition=='HATD1'"),
    x='Location',
    y='MAP',
    hue='Condition',
    hue_order=['LTD5', 'HATD1'],
    palette = ['C0','C4'],
    errorbar='se'
    )

#plt.ylim(.023,.026)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/anxietyReplay_MAP.pdf')

# %%
sns.lineplot(
    data=df.query("Condition=='HATDSwitch' or Condition=='LTD1'"),
    x='Location',
    y='MAP',
    palette=['C0','C4'],
    hue='Condition',
    errorbar='se'
    )
plt.plot([40,55],[35,35],'k')
plt.title("Switch day")
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/HATDS_MAP_location.pdf')
#%%

#%% SeqNMF
results_dir = '../../output_REM/seqNMF'
resultsList=os.listdir(results_dir)

data_list = []
for file_name in resultsList:
    if (
        file_name.startswith('seqReplayResults_') and file_name.endswith('.h5')
        and "pv1254" not in file_name
        ): # Exclude pv1254: # Only include seqResults, not replay results
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
df_replay_stats=df_replay.melt(id_vars=['state_ref','state_pred','mouse', 'condition'],value_name='numSeqs',value_vars=['S1_numSeqs', 'S2_numSeqs'],var_name='seqType')
df_replay_scores=df_replay.melt(id_vars=['state_ref','state_pred','mouse', 'condition'],value_name='seqScore',value_vars=['S1_score', 'S2_score'],var_name='seqType')

# %% Plot num. sequences vs anxiety
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_replay_stats.query("condition == 'LTD5' or condition == 'HATD1' and state_ref == 'wake' and state_pred == 'REMpost'"),
    hue='condition',
    y='numSeqs',
    x='seqType',
    palette=(['gray','C4']),
    # showfliers=False
    errorbar='se',
    capsize=.2
)

plt.xticks([0,1],['S1', 'S2'])
plt.ylabel('Num. significant \nsequences')
plt.xlabel('')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/anxiety_seqnmf_numSeqs_by_seqType.pdf')

#%% DESCRIPTIVES
#TODO

#%% STATS
#TODO
pg.rm_anova(data=df_replay_stats.query("seqType=='S1_numSeqs' and condition == 'LTD1' or condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'"),
         dv='numSeqs',
         within='condition',
         subject='mouse',
         )

# %% Plot seq. score vs anxiety
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_replay_scores.query("condition == 'LTD5' or condition == 'HATD1' and state_ref == 'wake' and state_pred == 'REMpost'"),
    x='condition',
    y='seqScore',
    hue='seqType',
    #palette=(['gray','C4']),
    # showfliers=False
    errorbar='se',
    capsize=.2
)

plt.xticks([0,1],['Control', 'Anxiety'], rotation=90)
plt.ylabel('Seq. score (z)')
plt.xlabel('')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/anxiety_seqnmf_seqScore_by_seqType.pdf')





# # %% Plot num. sequences vs novelty
# plt.figure(figsize=(.75,1))
# sns.barplot(
#     data=df_replay_scores.query("condition == 'LTD1' or condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'"),
#     x='condition',
#     y='seqScore',
#     palette=(['C3','gray']),
#     # showfliers=False
#     errorbar='se',
#     capsize=.2
# )

# #plt.xticks([0,1],['S1', 'S2'])
# plt.ylabel('Sequence score (z)')
# plt.xticks([0,1],['Novel','Familiar'])
# plt.xlabel('')
# plt.xticks(rotation=90)
# plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
# plt.savefig('../../output_REM/novely_seqnmf_seqScore.pdf')

# #%% STATS
# pg.rm_anova(data=df_replay_scores.query("condition == 'LTD1' or condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'"),
#          dv='seqScore',
#          within='condition',
#          subject='mouse',
#          )
# # %%

# %%
