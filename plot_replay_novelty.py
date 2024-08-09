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

#%% Plot examples
mouse = 'pv1060'
data_LTD1 = load_data(mouse='pv1060', condition='LTD1', state='REMpost', params=params)
data_LTD5 = load_data(mouse='pv1060', condition='LTD5', state='REMpost', params=params)

with h5py.File(os.path.join(params['path_to_output'],"neuron_selection", "selected_neurons_LTD1_pv1060.h5")) as f:
    selected_neurons_LTD1 = f['place_cells'][()]
with h5py.File(os.path.join(params['path_to_output'],"neuron_selection", "selected_neurons_LTD5_pv1060.h5")) as f:
    selected_neurons_LTD5 = f['place_cells'][()]

# Assembly stuff
with h5py.File(os.path.join(params['path_to_output'],"assembly", "assembly_LTD1_pv1060.h5")) as h5_file:
    LTD1_ID = h5_file['post_rem_A_ID'][()]
    LTD1_idx = h5_file['post_rem_react_idx'][()]
                
with h5py.File(os.path.join(params['path_to_output'],"assembly", "assembly_LTD5_pv1060.h5")) as h5_file:
    LTD5_ID = h5_file['post_rem_A_ID'][()]
    LTD5_idx = h5_file['post_rem_react_idx'][()]         

plt.figure(figsize=(4,1)) # TODO plot examples with Eric
plt.subplot(221)
plt.title('REM post (novel)')
plt.imshow(
     data_LTD1['binaryData'][:,selected_neurons_LTD1].T,
     aspect='auto',
     interpolation='none',
     cmap='gray_r',
     rasterized=True
     )
plt.xlim(0,3000)
plt.xticks([0,1800,3600],['','',''])
plt.yticks([0,128,256])
plt.ylabel('Neuron ID')

plt.subplot(425)
#event_colors = ['C{}'.format(i) for i in range(len(LTD1_ID))]
#plt.eventplot(np.stack((LTD1_idx,LTD1_ID),axis=1), colors=event_colors)
plt.eventplot(LTD1_idx[LTD1_ID==10], colors='C3')
plt.yticks([])
plt.xlim(0,3600)
plt.xticks([0,1800,3600],[0,1,2])
plt.xlabel('Time (min)')
plt.ylabel('')

plt.subplot(222)
plt.title('REM post (familiar)')
plt.imshow(
     data_LTD5['binaryData'][:,selected_neurons_LTD5].T,
     aspect='auto',
     interpolation='none',
     cmap='gray_r',
     rasterized=True
     )
plt.xlim(0,3600)
plt.xticks([0,1800,3600],['','',''])
plt.yticks([0,128,256],['','',''])
plt.ylabel('')

plt.subplot(426)
plt.eventplot(LTD5_idx[LTD5_ID==10], colors='gray')
plt.xlim(0,3600)
plt.xticks([0,1800,3600],[0,1,2])
plt.xlabel('Time (min)')
plt.yticks([])
plt.savefig("../../output_REM/replayNovelty_exampleAssemblies.pdf")

#%% Plot Bayesian decoding example
# Load posteriors
with h5py.File(os.path.join(params['path_to_output'], "posterior_probs", 'posterior_probs_LTD1_pv1060.h5'), 'r') as f:
    novel_posteriors = f[f'REMpost_posterior_probs'][()]

with h5py.File(os.path.join(params['path_to_output'], "posterior_probs", 'posterior_probs_LTD5_pv1060.h5'), 'r') as f:
    familiar_posteriors = f[f'REMpost_posterior_probs'][()]

# Load replay timestamps
with h5py.File(os.path.join(params['path_to_output'], "bayesian_replay", 'bayesian_replay_LTD1_pv1060_REMpost.h5'), 'r') as f:
    novel_replay_ts = f['replay_locs'][()]

with h5py.File(os.path.join(params['path_to_output'], "bayesian_replay", 'bayesian_replay_LTD5_pv1060_REMpost.h5'), 'r') as f:
    familiar_replay_ts = f['replay_locs'][()]

# Plot results
plt.figure(figsize=(4,1)) # TODO plot examples with Eric
plt.subplot(221)
plt.title('REM post (novel)')
plt.imshow(
     novel_posteriors.T,
     aspect='auto',
     interpolation='none',
     cmap='magma',
     rasterized=True,
     vmin=0,
     vmax=.075
     )
plt.xlim(0,3000)
plt.xticks([0,1800,3600],['','',''])
plt.yticks([0,19,39],[0,50,100])
plt.ylabel('Location (cm)')

plt.subplot(425)
plt.eventplot(novel_replay_ts,colors='C3')
plt.yticks([])
plt.xlim(0,3600)
plt.xticks([0,1800,3600],[0,1,2])
plt.xlabel('Time (min)')
plt.ylabel('')

plt.subplot(222)
plt.title('REM post (familiar)')
plt.imshow(
     familiar_posteriors.T,
     aspect='auto',
     interpolation='none',
     cmap='magma',
     rasterized=True,
     vmin=0,
     vmax=.075
     )
plt.xlim(0,3600)
plt.xticks([0,1800,3600],['','',''])
plt.yticks([0,19,39],['','',''])
plt.ylabel('')

plt.subplot(426)
plt.eventplot(familiar_replay_ts,colors='gray')
plt.xlim(0,3600)
plt.xticks([0,1800,3600],[0,1,2])
plt.xlabel('Time (min)')
plt.yticks([])
plt.savefig("../../output_REM/replayNovelty_exampleBayesianReplay.pdf")


#%% Plot seqNMF examples
# load replay timestamps
with h5py.File(
                os.path.join(
                    params["path_to_output"],
                    'seqNMF_timestamps',
                    'seqReplayResults_LTD1_pv1060_wake_REMpost.h5'
                ),
                "r",
            ) as seqNMF_file:
                LTD1_seqNMF_ts = seqNMF_file['seqReplayLocs'][()]
                
with h5py.File(
                os.path.join(
                    params["path_to_output"],
                    'seqNMF_timestamps',
                    'seqReplayResults_LTD5_pv1060_wake_REMpost.h5'
                ),
                "r",
            ) as seqNMF_file:
                LTD5_seqNMF_ts = seqNMF_file['seqReplayLocs'][()]

plt.figure(figsize=(4,1)) # TODO plot examples with Eric
plt.subplot(221)
plt.title('REM post (novel)')
plt.imshow(
     data_LTD1['binaryData'][:,selected_neurons_LTD1].T,
     aspect='auto',
     interpolation='none',
     cmap='gray_r',
     rasterized=True
     )
plt.xlim(0,3000)
plt.xticks([0,1800,3600],['','',''])
plt.yticks([0,128,256])
plt.ylabel('Neuron ID')

plt.subplot(425)
#event_colors = ['C{}'.format(i) for i in range(len(LTD1_ID))]
#plt.eventplot(np.stack((LTD1_idx,LTD1_ID),axis=1), colors=event_colors)
plt.eventplot(LTD1_seqNMF_ts, colors='C3')
plt.yticks([])
plt.xlim(0,3600)
plt.xticks([0,1800,3600],[0,1,2])
plt.xlabel('Time (min)')
plt.ylabel('')

plt.subplot(222)
plt.title('REM post (familiar)')
plt.imshow(
     data_LTD5['binaryData'][:,selected_neurons_LTD5].T,
     aspect='auto',
     interpolation='none',
     cmap='gray_r',
     rasterized=True
     )
plt.xlim(0,3600)
plt.xticks([0,1800,3600],['','',''])
plt.yticks([0,128,256],['','',''])
plt.ylabel('')

plt.subplot(426)
plt.eventplot(LTD5_seqNMF_ts, colors='gray')
plt.xlim(0,3600)
plt.xticks([0,1800,3600],[0,1,2])
plt.xlabel('Time (min)')
plt.yticks([])
plt.savefig("../../output_REM/replayNovelty_exampleSeqNMF.pdf")


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
        and "pv1252" not in file_name # Exclude pv1254
    ):
        h5_file = h5py.File(os.path.join(results_dir, file_name))
        data = load_data(h5_file["mouse"][0].decode("utf-8"),
                         h5_file["condition"][0].decode("utf-8"),
                         'REMpost',
                         params
                         )
        numFrames = data['binaryData'].shape[0]
        recordingLength = numFrames/params['sampling_frequency']

        # Percent replay?
        ID_wake_assemblies = np.unique(h5_file['post_rem_A_ID'][()])
        ID_replayed_assemblies = np.unique(h5_file['post_rem_A_ID'][()][h5_file['post_rem_A_sig'][()]==1])
        num_wake_assemblies = len(np.unique(h5_file["post_rem_A_ID"][()]))
        num_REMpost_assemblies = len(np.unique(h5_file['post_rem_A_ID'][()][h5_file['post_rem_A_sig'][()]==1]))
        portion_replayed = num_REMpost_assemblies/num_wake_assemblies

        data_list_num_events.append(
            {
                "mouse": h5_file["mouse"][0].decode("utf-8"),
                "condition": h5_file["condition"][0].decode("utf-8"),
                "numSigEvents": np.sum(h5_file['post_rem_A_sig'][()]==1),
                "numAssemblies": np.max(h5_file['post_rem_A_ID']),
                "meanReactPerAssembly": np.mean(np.unique(h5_file['post_rem_A_ID'][()],return_counts=True)),
                "meanFreqPerAssembly": np.mean(np.unique(h5_file['post_rem_A_ID'][()],return_counts=True))/recordingLength,
                "numWakeAssemblies": num_wake_assemblies,
                "numREMpostAssemblies": num_REMpost_assemblies,
                "portionREMpostReplay": portion_replayed
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
                }
            )

        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)
df_numEvents = pd.DataFrame(data_list_num_events)    

#%% Plot awake assemblies
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_numEvents.query("condition=='LTD1' or condition=='LTD5'"),
    x='condition',
    y='numWakeAssemblies',
    order=['LTD1', 'LTD5'],
    palette=(['C3','gray']),
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_numEvents.query("condition=='LTD1' or condition=='LTD5'"),
    color='k',
    order=['LTD1', 'LTD5'],
    x='condition',
    y='numWakeAssemblies',
    size=2
)
plt.xticks([0,1],['Novel','Familiar'],rotation=90)
plt.xlabel('')
plt.ylabel('Portion replayed\nassemblies')
plt.savefig("../../output_REM/novelty_numAwakeAssemblies.pdf")

#%% Plot percent replayed
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_numEvents.query("condition=='LTD1' or condition=='LTD5'"),
    x='condition',
    y='portionREMpostReplay',
    order=['LTD1', 'LTD5'],
    palette=(['C3','gray']),
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_numEvents.query("condition=='LTD1' or condition=='LTD5'"),
    color='k',
    order=['LTD1', 'LTD5'],
    x='condition',
    y='portionREMpostReplay',
    size=2
)
plt.xticks([0,1],['Novel','Familiar'],rotation=90)
plt.xlabel('')
plt.ylabel('Portion replayed\nassemblies')
plt.savefig("../../output_REM/novelty_portionReplayedAssemblies.pdf")

#%% Plot results
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_numEvents.query("condition=='LTD1' or condition=='LTD5'"),
    x='condition',
    y='numAssemblies',
    order=['LTD1', 'LTD5'],
    palette=(['C3','gray']),
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_numEvents.query("condition=='LTD1' or condition=='LTD5'"),
    color='k',
    order=['LTD1', 'LTD5'],
    x='condition',
    y='numAssemblies',
    size=2
)
plt.xticks([0,1],['Novel','Familiar'],rotation=90)
plt.xlabel('')
plt.ylabel('Num. replayed\nassemblies')
plt.savefig("../../output_REM/noveltyReplay_numAssemblies.pdf")

#%% DESCRIPTIVES
mean_novel = df_numEvents.query("condition=='LTD1'")['numAssemblies'].mean()
SEM_novel = df_numEvents.query("condition=='LTD1'")['numAssemblies'].sem()
print(f'{mean_novel} +/- {SEM_novel}')

mean_fam = df_numEvents.query("condition=='LTD5'")['numAssemblies'].mean()
SEM_fam = df_numEvents.query("condition=='LTD5'")['numAssemblies'].sem()
print(f'{mean_fam} +/- {SEM_fam}')

#%% STATS
pg.ttest(
    x=df_numEvents.query("condition=='LTD1'")['numAssemblies'],
    y=df_numEvents.query("condition=='LTD5'")['numAssemblies'],
    paired=True)

#%% Plot frequency
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_numEvents.query("condition=='LTD1' or condition=='LTD5'"),
    x='condition',
    y='meanFreqPerAssembly',
    order=['LTD1', 'LTD5'],
    palette=(['C3','gray']),
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_numEvents.query("condition=='LTD1' or condition=='LTD5'"),
    color='k',
    order=['LTD1', 'LTD5'],
    x='condition',
    y='meanFreqPerAssembly',
    size=2
)

plt.xticks([0,1],['Novel','Familiar'],rotation=90)
plt.xlabel('')
plt.ylabel('Mean reactivation\nfrequency (Hz)')
plt.savefig("../../output_REM/noveltyReplay_meanAssemblyFreq.pdf")

#%% DESCRIPTIVES
mean_novel = df_numEvents.query("condition=='LTD1'")['meanFreqPerAssembly'].mean()
SEM_novel = df_numEvents.query("condition=='LTD1'")['meanFreqPerAssembly'].sem()
print(f'{mean_novel} +/- {SEM_novel}')

mean_fam = df_numEvents.query("condition=='LTD5'")['meanFreqPerAssembly'].mean()
SEM_fam = df_numEvents.query("condition=='LTD5'")['meanFreqPerAssembly'].sem()
print(f'{mean_fam} +/- {SEM_fam}')

#%% STATS
pg.ttest(
    x=df_numEvents.query("condition=='LTD1'")['meanFreqPerAssembly'],
    y=df_numEvents.query("condition=='LTD5'")['meanFreqPerAssembly'],
    paired=True)

#%%
# plt.figure(figsize=(1.25,1))
# sns.boxenplot(
#     data=df.query("condition=='LTD1' or condition=='LTD5'"),
#     x = 'condition',
#     y ='replayEventStrength',
#     order=['LTD1', 'LTD5'],
#     palette=(['C3','gray']),
#     showfliers=False
# )
# # plt.title('REM post')
# plt.xlabel('')
# plt.ylabel('Strength (z)')
# plt.yscale('log')
# plt.savefig("../../output_REM/noveltyReplay_assemblyStrength.pdf")

#%% Bayesian replay
results_dir = params['path_to_output']+"/bayesian_replay"
resultsList = os.listdir(results_dir)

data_list = []
num_event_list = []

for file_name in tqdm(resultsList):
    if (
        file_name.startswith("bayesian_replay_")
        and file_name.endswith(".h5")
        and "pv1252" not in file_name
    ):
        h5_file = h5py.File(os.path.join(results_dir, file_name))
        data = load_data(h5_file["mouse"][()].decode("utf-8"),
                         h5_file["condition"][()].decode("utf-8"),
                         'REMpost',
                         params
                         )
        numFrames = data['binaryData'].shape[0]
        recordingLength = numFrames/params['sampling_frequency']

        num_event_list.append(
            {
                "mouse": h5_file["mouse"][()].decode("utf-8"),
                "condition": h5_file["condition"][()].decode("utf-8"),
                "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                "numReplayEvents": len(h5_file["replay_locs"][()]),
                "eventsFrequency": len(h5_file["replay_locs"][()])/recordingLength
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
df_numEvents = df_numEvents.query("Type=='replay'")

#%% Plot number of events per condition
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_numEvents.query("condition=='LTD1' or condition=='LTD5'"),
    x='condition',
    y='eventsFrequency',
    order=['LTD1', 'LTD5'],
    palette=(['C3','gray']),
    errorbar='se',
    capsize=.2
)

sns.stripplot(
    data=df_numEvents.query("condition=='LTD1' or condition=='LTD5'"),
    x='condition',
    y='eventsFrequency',
    order=['LTD1', 'LTD5'],
    color='k',
    size=2
)

plt.ylabel('Replay events\nfrequency (Hz)')
plt.xticks([0,1],['Novel','Familiar'])
plt.xlabel('')
plt.xticks(rotation=90)
#plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/novelty_bayesian_eventsFrequency.pdf')

#%% DESCRIPTIVES
mean_novelty = df_numEvents.query("condition=='LTD1'")['eventsFrequency'].mean()
SEM_novelty = df_numEvents.query("condition=='LTD1'")['eventsFrequency'].sem()
print(f'Novel: {mean_novelty} +/- {SEM_novelty}')

mean_familiar = df_numEvents.query("condition=='LTD5'")['eventsFrequency'].mean()
SEM_familiar = df_numEvents.query("condition=='LTD5'")['eventsFrequency'].sem()
print(f'Familiar: {mean_familiar} +/- {SEM_familiar}')

#%% STATS
pg.ttest(
    x=df_numEvents.query("condition=='LTD1'")['eventsFrequency'],
    y=df_numEvents.query("condition=='LTD5'")['eventsFrequency'],
    paired=True
)

#%% Plot replay score as a function of novelty
plt.figure(figsize=(.75,1))
sns.boxenplot(
    data=df.query("condition=='LTD1' or condition=='LTD5'"),
    x='condition',
    y='replayEventScore',
    order=['LTD1', 'LTD5'],
    palette=(['C3','gray']),
    showfliers=False
    #units='mouse',
    #errorbar='se',
    #capsize=.2
)

plt.ylabel('Replay score ($R^{2}$)')
plt.xticks([0,1],['Novel','Familiar'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/novelty_bayesian_replay_scores.pdf')

#%% DESCRIPTIVES
mean_novelty = df.query("condition=='LTD1'")['replayEventScore'].mean()
SEM_novelty = df.query("condition=='LTD1'")['replayEventScore'].sem()
print(f'Novel: {mean_novelty} +/- {SEM_novelty}')

mean_familiar = df.query("condition=='LTD5'")['replayEventScore'].mean()
SEM_familiar = df.query("condition=='LTD5'")['replayEventScore'].sem()
print(f'Familiar: {mean_familiar} +/- {SEM_familiar}')

#%% STATS
pg.ttest(
    x=df.query("condition=='LTD1'")['replayEventScore'],
    y=df.query("condition=='LTD5'")['replayEventScore'],
    paired=False
)

#%% Plot jumpiness as a function of novelty
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df.query("condition=='LTD1' or condition=='LTD5'"),
    x='condition',
    y='replayEventJumpiness',
    order=['LTD1', 'LTD5'],
    palette=(['C3','gray']),
    #showfliers=False
    errorbar='se',
    capsize=.2
)

plt.ylabel('Max. jumpiness (cm)')
plt.xticks([0,1],['Novel','Familiar'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/novelty_bayesian_replay_jumpiness.pdf')

#%% DESCRIPTIVES
mean_novelty = df.query("Type=='replay' and condition=='LTD1'")['replayEventJumpiness'].mean()
SEM_novelty = df.query("Type=='replay' and condition=='LTD1'")['replayEventJumpiness'].sem()
print(f'Novel: {mean_novelty} +/- {SEM_novelty}')

mean_familiar = df.query("Type=='replay' and condition=='LTD5'")['replayEventJumpiness'].mean()
SEM_familiar = df.query("Type=='replay' and condition=='LTD5'")['replayEventJumpiness'].sem()
print(f'Familiar: {mean_familiar} +/- {SEM_familiar}')

#%% STATS
pg.anova(
    data=df.query("Type=='replay' and condition=='LTD1' or condition=='LTD5'"),
    dv='replayEventJumpiness',
    between='condition',
)

#%%
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df.query("Type=='replay' and condition=='LTD1' or condition=='LTD5'"),
    x='condition',
    y='replayEventSlope',
    order=['LTD1', 'LTD5'],
    palette=(['C3','gray']),
    # showfliers=False
    errorbar='se',
    capsize=.2
)

plt.ylabel('Slope (cm.s$^{-1}$)')
plt.xticks([0,1],['Novel','Familiar'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/novelty_bayesian_replay_slope.pdf')


#%% SeqNMF
results_dir = '../../output_REM/seqNMF'
resultsList=os.listdir(results_dir)

#%%
data_list = []
for file_name in resultsList:
    if (
        file_name.startswith('seqReplayResults_')
        and file_name.endswith('.h5')
        and "pv1252" not in file_name
        ):
        h5_file = h5py.File(os.path.join(results_dir,file_name), 'r')
        data = load_data(h5_file["mouse"][()].decode("utf-8"),
                         h5_file["condition"][()].decode("utf-8"),
                         'REMpost',
                         params
                         )
        numFrames = data['binaryData'].shape[0]
        recordingLength = numFrames/params['sampling_frequency']

        # Per event analysis
        # for i in range(len(h5_file["replay_locs"][()])):
        #     data_list.append(  # This will create one list entry per cell
        #         {
        #             "eventID": i,
        #             "mouse": h5_file["mouse"][()].decode("utf-8"),
        #             "condition": h5_file["condition"][()].decode("utf-8"),
        #             "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
        #             "replayEventTime": h5_file['replay_locs'][i]/params['sampling_frequency'],
        #             "replayEventScore": h5_file['replay_score'][i],
        #             "replayEventJumpiness": h5_file['replay_jumpiness'][i],
        #             "replayEventLength": h5_file['replay_length'][i],
        #             "replayEventSlope": abs(h5_file['replay_slope'][i])
        #         }
        #     )




        data_list.append( #This will create one list entry per cell
                {
                    'state_ref':h5_file['state_ref'][()].decode("utf-8"),
                    'state_pred':h5_file['state_pred'][()].decode("utf-8"),
                    'condition':h5_file['condition'][()].decode("utf-8"),
                    'mouse':h5_file['mouse'][()].decode("utf-8"),
                    'S1_numSeqs':h5_file['S1_numSeqs'][()],
                    'S2_numSeqs':h5_file['S2_numSeqs'][()],
                    'S1_freq':h5_file['S1_numSeqs'][()]/recordingLength,
                    'S2_freq':h5_file['S2_numSeqs'][()]/recordingLength,
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
df_replay_freq=df_replay.melt(id_vars=['state_ref','state_pred','mouse', 'condition'],value_name='seqFreq',value_vars=['S1_freq', 'S2_freq'],var_name='seqType')
df_replay_scores=df_replay.melt(id_vars=['state_ref','state_pred','mouse', 'condition'],value_name='seqScore',value_vars=['S1_score', 'S2_score'],var_name='seqType')

# %% Plot num. sequences vs novelty
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_replay_stats.query("seqType=='S1_numSeqs'").query("condition == 'LTD1' or condition == 'LTD5'").query("state_ref == 'wake' and state_pred == 'REMpost'"),
    x='condition',
    y='numSeqs',
    palette=(['C3','gray']),
    #alpha=.2,
    #linewidth=0
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_replay_stats.query("seqType=='S1_numSeqs'").query("condition == 'LTD1' or condition == 'LTD5'").query("state_ref == 'wake' and state_pred == 'REMpost'"),
    x='condition',
    y='numSeqs',
    palette=(['k','k']),
    size=2
    # errorbar='se',
    # capsize=.2
)

#plt.xticks([0,1],['S1', 'S2'])
plt.ylabel('Num. significant \nsequences')
plt.xticks([0,1],['Novel','Familiar'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/novely_seqnmf_numSeqs.pdf')

#%% DESCRIPTIVES
mean_novelty = df_replay_stats.query("seqType=='S1_numSeqs'").query("condition == 'LTD1' and state_ref == 'wake' and state_pred == 'REMpost'")['numSeqs'].mean()
SEM_novelty = df_replay_stats.query("seqType=='S1_numSeqs' and condition == 'LTD1' and state_ref == 'wake' and state_pred == 'REMpost'")['numSeqs'].sem()
print(f'Novel: {mean_novelty} +/- {SEM_novelty}')

mean_familiar = df_replay_stats.query("seqType=='S1_numSeqs'").query("condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'")['numSeqs'].mean()
SEM_familiar = df_replay_stats.query("seqType=='S1_numSeqs' and condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'")['numSeqs'].sem()
print(f'Familiar: {mean_familiar} +/- {SEM_familiar}')

#%% STATS
pg.rm_anova(data=df_replay_stats.query("seqType=='S1_numSeqs'").query("condition == 'LTD1' or condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'"),
         dv='numSeqs',
         within='condition',
         subject='mouse',
         )

#%% DESCRIPTIVES
mean_novelty = df_replay_freq.query("seqType=='S1_freq' and condition == 'LTD1' and state_ref == 'wake' and state_pred == 'REMpost'")['seqFreq'].mean()
SEM_novelty = df_replay_freq.query("seqType=='S1_freq' and condition == 'LTD1' and state_ref == 'wake' and state_pred == 'REMpost'")['seqFreq'].sem()
print(f'Novel: {mean_novelty} +/- {SEM_novelty}')

mean_familiar = df_replay_freq.query("seqType=='S1_freq' and condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'")['seqFreq'].mean()
SEM_familiar = df_replay_freq.query("seqType=='S1_freq' and condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'")['seqFreq'].sem()
print(f'Familiar: {mean_familiar} +/- {SEM_familiar}')

#%% STATS
pg.rm_anova(data=df_replay_freq.query("seqType=='S1_freq' and condition == 'LTD1' or condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'"),
         dv='seqFreq',
         within='condition',
         subject='mouse',
         )


#%%
results_dir = '../../output_REM/seqNMF_timestamps'
resultsList=os.listdir(results_dir)

#%% Focus on replayu timestamps to compute instantaneous frequencies
data_list = []
for file_name in resultsList:
    if (
        file_name.startswith('seqReplayResults_')
        and file_name.endswith('.h5')
        and "pv1252" not in file_name
        ):
        h5_file = h5py.File(os.path.join(results_dir,file_name), 'r')
        data = load_data(h5_file["mouse"][()].decode("utf-8"),
                         h5_file["condition"][()].decode("utf-8"),
                         'REMpost',
                         params
                         )
        numFrames = data['binaryData'].shape[0]
        recordingLength = numFrames/params['sampling_frequency']
        inst_frequencies = h5_file['seqReplayLocs'][()]/params['sampling_frequency']
        inst_frequencies = 1/np.diff(inst_frequencies)

        # Per event analysis
        for i in range(len(inst_frequencies)):
            data_list.append(  # This will create one list entry per cell
                {
                    "eventID": i,
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    'state_ref':h5_file['state_ref'][()].decode("utf-8"),
                    'state_pred':h5_file['state_pred'][()].decode("utf-8"),
                    "inst_freq": inst_frequencies[i]
                }
            )

        # Close files
        h5_file.close()

df_replay = pd.DataFrame(data_list)

#%%
plt.figure(figsize=(.75,1))
sns.violinplot(
    data=df_replay.query("condition == 'LTD1' or condition == 'LTD5'").query("state_ref == 'wake' and state_pred == 'REMpost'"),
    x='condition',
    y='inst_freq',
    palette=(['C3','gray']),
    linewidth=0,
    alpha=.2
    # errorbar='se',
    # capsize=.2
)
sns.stripplot(
    data=df_replay.query("condition == 'LTD1' or condition == 'LTD5'").query("state_ref == 'wake' and state_pred == 'REMpost'"),
    x='condition',
    y='inst_freq',
    palette=(['C3','gray']),
    size=2
    # errorbar='se',
    # capsize=.2
)

#plt.xticks([0,1],['S1', 'S2'])
plt.ylabel('Num. significant \nsequences')
plt.xticks([0,1],['Novel','Familiar'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/novely_seqnmf_inst_freq.pdf')







# %% Plot num. sequences vs novelty
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_replay_scores.query("condition == 'LTD1' or condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'"),
    x='condition',
    y='seqScore',
    palette=(['C3','gray']),
    # showfliers=False
    errorbar='se',
    capsize=.2
)

#plt.xticks([0,1],['S1', 'S2'])
plt.ylabel('Sequence score (z)')
plt.xticks([0,1],['Novel','Familiar'])
plt.xlabel('')
plt.xticks(rotation=90)
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/novely_seqnmf_seqScore.pdf')

#%% STATS
pg.rm_anova(data=df_replay_scores.query("condition == 'LTD1' or condition == 'LTD5' and state_ref == 'wake' and state_pred == 'REMpost'"),
         dv='seqScore',
         within='condition',
         subject='mouse',
         )
# %%    
