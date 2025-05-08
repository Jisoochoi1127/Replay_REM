# %% Imports
import numpy as np
import os
import pandas as pd
import pingouin as pg
import seaborn as sns
import matplotlib.pyplot as plt
import h5py
import yaml
from tqdm import tqdm
from utils.helperFunctions import load_data, extract_seq_score, extract_H

plt.style.use("plot_style.mplstyle")

# %% Load parameters
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

np.random.seed(params["seed"])

#%% Assembly analysis
results_dir = params['path_to_output']+"/assembly"
resultsList = os.listdir(results_dir)
data_list = []
data_list_num_events = []

for file_name in tqdm(resultsList):
    if (
        file_name.startswith("assembly_")
        and file_name.endswith(".h5")
        # and "pv1254" not in file_name # Exclude pv1254
    ):
        h5_file = h5py.File(os.path.join(results_dir, file_name))
        data = load_data(h5_file["mouse"][0].decode("utf-8"),
                         h5_file["condition"][0].decode("utf-8"),
                         'REMpre',
                         params
                         )
        numFrames = data['binaryData'].shape[0]
        recordingLength = numFrames/params['sampling_frequency']

        data_list_num_events.append(
            {
                "mouse": h5_file["mouse"][0].decode("utf-8"),
                "condition": h5_file["condition"][0].decode("utf-8"),
                "numSigEvents": np.sum(h5_file['pre_rem_A_sig'][()]==1),
                "numPreplayedAssemblies": np.max(h5_file['pre_rem_A_ID']),
                "numReplayedAssemblies": np.max(h5_file['post_rem_A_ID']),
                "meanReactPerAssembly": np.mean(np.unique(h5_file['pre_rem_A_ID'][()],return_counts=True)),
                "meanPreplayFreqPerAssembly": np.mean(np.unique(h5_file['pre_rem_A_ID'][()],return_counts=True))/recordingLength,
                "meanReplayFreqPerAssembly": np.mean(np.unique(h5_file['post_rem_A_ID'][()],return_counts=True))/recordingLength
            }
        )

        for i in range(len(h5_file["pre_rem_A_react_idx"][()])):
            data_list.append(  # This will create one list entry per cell
                {
                    "eventID": i,
                    "mouse": h5_file["mouse"][0].decode("utf-8"),
                    "condition": h5_file["condition"][0].decode("utf-8"),
                    "replayEventTime": h5_file['pre_rem_A_react_idx'][i]/params['sampling_frequency'],
                    "replayEventID": h5_file['pre_rem_A_ID'][i],
                    "replayEventStrength": h5_file['pre_rem_A_str'][i],
                    "replayEventSignificance": h5_file['pre_rem_A_sig'][i],
                }
            )

        # Close files
        h5_file.close()

df = pd.DataFrame(data_list)
df_numEvents = pd.DataFrame(data_list_num_events)

#%% Melt data
melted_df_numEvents = df_numEvents.melt(
    id_vars = ['mouse', 'condition'],
    value_vars = ['numPreplayedAssemblies', 'numReplayedAssemblies'],
    var_name = 'Type',
    value_name = 'Num. assemblies'
)

#%% Plot results
plt.figure(figsize=(1.25,.5))
sns.barplot(
    data=melted_df_numEvents.query("condition=='LTD1'"),
    y='Type',
    hue='Type',
    palette=['C3','C0'],
    x='Num. assemblies',
    errorbar='se',
    capsize=.2,
    legend=False
)
sns.stripplot(
    data=melted_df_numEvents.query("condition=='LTD1'"),
    hue='Type',
    y='Type',
    palette=['gray','gray'],
    x='Num. assemblies',
    size=2,
    legend=False
)
plt.xlim(0,25)
plt.yticks([0,1],['Pre-task REM','Post-task REM'])
plt.ylabel('')
#plt.xlabel('Number of\npre-existing assemblies')
#plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)

plt.savefig("../../output_REM/REMpre_numAssemblies.pdf")

#%% DESCRIPTIVES
mean_preplay = df_numEvents.query("condition=='LTD1'")['numPreplayedAssemblies'].mean()
SEM_preplay = df_numEvents.query("condition=='LTD1'")['numPreplayedAssemblies'].sem()
print(f'{mean_preplay} +/- {SEM_preplay}')

mean_replayed = df_numEvents.query("condition=='LTD1'")['numReplayedAssemblies'].mean()
SEM_replayed = df_numEvents.query("condition=='LTD1'")['numReplayedAssemblies'].sem()
print(f'{mean_replayed} +/- {SEM_replayed}')

#%% STATS
pg.ttest(
    x=df_numEvents.query("condition=='LTD1'")['numPreplayedAssemblies'],
    y=df_numEvents.query("condition=='LTD1'")['numReplayedAssemblies'],
    paired=True)

#%% Assembly reactivations frequency
melted_df_freq = df_numEvents.melt(
    id_vars = ['mouse', 'condition'],
    value_vars = ['meanPreplayFreqPerAssembly', 'meanReplayFreqPerAssembly'],
    var_name = 'Type',
    value_name = 'Reactivation frequency (Hz)'
)

#%%
plt.figure(figsize=(1.25,.5))
sns.barplot(
    data=melted_df_freq.query("condition=='LTD1'"),
    y='Type',
    hue='Type',
    palette=['C3','C0'],
    x='Reactivation frequency (Hz)',
    errorbar='se',
    capsize=.2,
    legend=False
)
sns.stripplot(
    data=melted_df_freq.query("condition=='LTD1'"),
    hue='Type',
    y='Type',
    palette=['gray','gray'],
    x='Reactivation frequency (Hz)',
    size=2,
    legend=False
)
#plt.xlim(0,25)
plt.yticks([0,1],['Pre-task REM','Post-task REM'])
plt.xlabel('Reactivation\nfrequency (Hz)')
plt.ylabel('')
#plt.xlabel('Number of\npre-existing assemblies')
#plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)

plt.savefig("../../output_REM/REMpre_reac_freq.pdf")

#%% Descriptives
mean_preplay = df_numEvents.query("condition=='LTD1'")['meanPreplayFreqPerAssembly'].mean()
SEM_preplay = df_numEvents.query("condition=='LTD1'")['meanPreplayFreqPerAssembly'].sem()
print(f'{mean_preplay} +/- {SEM_preplay}')

mean_replayed = df_numEvents.query("condition=='LTD1'")['meanReplayFreqPerAssembly'].mean()
SEM_replayed = df_numEvents.query("condition=='LTD1'")['meanReplayFreqPerAssembly'].sem()
print(f'{mean_replayed} +/- {SEM_replayed}')

#%% STATS
pg.ttest(
    x=df_numEvents.query("condition=='LTD1'")['meanPreplayFreqPerAssembly'],
    y=df_numEvents.query("condition=='LTD1'")['meanReplayFreqPerAssembly'],
    paired=True)

#%% Load example
with h5py.File(os.path.join(params['path_to_output'],"neuron_selection", "selected_neurons_LTD1_pv1060.h5")) as f:
    selected_neurons = f['place_cells'][()]
# data_LTD1 = load_data(mouse='pv1060', condition='LTD1', state='wake', params=params)
data_REMpre = load_data(mouse='pv1060', condition='LTD1', state='REMpre', params=params)

#%% Plot example
plt.figure(figsize=(3,1))
plt.imshow(
    data_REMpre['binaryData'][:,selected_neurons].T,
    interpolation='none',
    aspect='auto',
    cmap='Blues',
    rasterized=True,
    vmin=0,
    vmax=1
)
event_times = df.query("mouse=='pv1060' and replayEventID==8 and condition=='LTD1'")['replayEventTime']

for key in event_times.keys():
    timeStamp = event_times[key]*30
    plt.plot([timeStamp, timeStamp],[-10,0],
             color='C4')
plt.xlim(0,2700)
plt.xticks([0,900,1800,2700],[0,30,60,90])
plt.xlabel('Time (s)')
plt.ylabel('Neuron ID')
plt.yticks([0,128,256])
plt.title('Pre-task REM (naive)')
plt.savefig("../../output_REM/preplay_assembly_example.pdf")

#%% SeqNMF
# Load example
mouse = 'pv1060'
with h5py.File(os.path.join(params['path_to_output'],"neuron_selection", "selected_neurons_LTD1_pv1060.h5")) as f:
    selected_neurons = f['place_cells'][()]
# data_LTD1 = load_data(mouse='pv1060', condition='LTD1', state='wake', params=params)
data_REMpre = load_data(mouse='pv1060', condition='LTD1', state='REMpre', params=params)

seqReplay_scores, seqReplay_pvalues, seqReplay_locs, W_ref = extract_seq_score(data_REMpre['binaryData'][:,selected_neurons],
                                                                                params)
#%% Sorting index
seqSortingIndex = np.argsort(np.argmax(np.concatenate((W_ref[:,0,:],W_ref[:,1,:]),axis=1),axis=1))

H = extract_H(W_ref, data_REMpre['binaryData'][:,selected_neurons])
SF1 = data_REMpre['binaryData'][:,selected_neurons]*H[0][:,None]
SF2 = data_REMpre['binaryData'][:,selected_neurons]*H[1][:,None]
sequence_data = SF1+SF2

#%%
plt.figure(figsize=(3,1))
plt.imshow(sequence_data[:,seqSortingIndex].T,
           interpolation='none',
           aspect='auto',
           cmap='Blues',
           rasterized=True,
           vmin=0,
           vmax=1
           )
plt.xlim(500,2700)
plt.xticks([500,1400,2300],[0,30,60])
plt.xlabel('Time (s)')
plt.ylabel('Neuron ID')
plt.yticks([0,128,256])
plt.title('Pre-task REM (naive)')
plt.savefig("../../output_REM/preplay_seqNMF_example.pdf")

#%% Analysis
results_dir = '../../output_REM/seqNMF'
resultsList=os.listdir(results_dir)

#%%
data_list = []
event_data_list = []
for file_name in resultsList:
    if file_name.startswith('seqReplayResults_') and file_name.endswith('.h5'): # Only include seqResults, not replay results
        h5_file = h5py.File(os.path.join(results_dir,file_name))
        data = load_data(h5_file["mouse"][()].decode("utf-8"),
                         h5_file["condition"][()].decode("utf-8"),
                         'REMpost',
                         params
                         )
        numFrames = data['binaryData'].shape[0]
        recordingLength = numFrames/params['sampling_frequency']
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
df_replay_stats=df_replay.melt(id_vars=['state_ref','state_pred','mouse', 'condition'],value_name='numSeqs',value_vars=['S1_numSeqs', 'S2_numSeqs'],var_name='seqType')

#%%
print("Avg seqNMF replay events: " + str(df_replay.query('condition=="LTD1"')['S1_freq'].mean()))
print("Avg seqNMF replay events: " + str(df_replay.query('condition=="LTD1"')['S1_freq'].sem()))
#%%
plt.figure(figsize=(1,.75))
sns.heatmap(data=df_replay_stats.query("condition=='LTD1'").pivot_table(
    index='state_ref',
    columns='state_pred', 
    values='numSeqs'
    ),
    cmap='magma',
    vmin=0,
    rasterized=True,
    cbar_kws={'label':'Num. replay events'}
)
#plt.xticks([])
#plt.yticks([])
plt.xlabel('Reference')
plt.ylabel('Target')
plt.savefig('../../output_REM/intStruct_numSeqs_heatmap.pdf')

#%% Plot seqNMF preplay
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_replay.query("condition=='LTD1' and state_ref=='REMpre' and state_pred=='wake'"),
    x=1,
    y='S1_numSeqs',
    # showfliers=False,
    color='C3',
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_replay.query("condition=='LTD1' and state_ref=='REMpre' and state_pred=='wake'"),
    x=1,
    y='S1_numSeqs',
    color='k',
    size=2,
    dodge=False,
    legend=False
)

sns.barplot(
    data=df_replay.query("condition=='LTD1' and state_ref=='wake' and state_pred=='REMpost'"),
    x=2,
    y='S1_numSeqs',
    # showfliers=False,
    color='C0',
    errorbar='se',
    capsize=.2
)
sns.stripplot(
    data=df_replay.query("condition=='LTD1' and state_ref=='wake' and state_pred=='REMpost'"),
    x=2,
    y='S1_numSeqs',
    color='k',
    size=2,
    dodge=False,
    legend=False
)

plt.xticks([0,1],['Preplay', 'Replay'],rotation=90)
plt.ylabel('Num. rigid\nsequences')
plt.xlabel('')
plt.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
plt.savefig('../../output_REM/seqNMF_preplay_numSeqs.pdf')

#%% DESCRIPTIVES
mean_preplay = df_replay.query("condition=='LTD1' and state_ref=='REMpre' and state_pred=='wake'")['S1_numSeqs'].mean()
SEM_preplay = df_replay.query("condition=='LTD1' and state_ref=='REMpre' and state_pred=='wake'")['S1_numSeqs'].sem()
print(f'{mean_preplay} +/- {SEM_preplay}')

mean_replay = df_replay.query("condition=='LTD1' and state_ref=='wake' and state_pred=='REMpost'")['S1_numSeqs'].mean()
SEM_replay = df_replay.query("condition=='LTD1' and state_ref=='wake' and state_pred=='REMpost'")['S1_numSeqs'].sem()
print(f'{mean_replay} +/- {SEM_replay}')

#%% STATS
pg.ttest(
    x=df_replay.query("condition=='LTD1' and state_ref=='REMpre' and state_pred=='wake'")['S1_numSeqs'],
    y=df_replay.query("condition=='LTD1' and state_ref=='wake' and state_pred=='REMpost'")['S1_numSeqs'],
    paired=True)

#%% Flexible replay
# Example
with h5py.File(os.path.join(params['path_to_output'], "posterior_probs", 'posterior_probs_LTD1_pv1060.h5'), 'r') as f:
    REMpre_posteriors = f[f'REMpre_posterior_probs'][()]

# Load preplay timestamps
with h5py.File(os.path.join(params['path_to_output'], "bayesian_replay", 'bayesian_replay_LTD1_pv1060_REMpre.h5'), 'r') as f:
    novel_replay_ts = f['replay_locs'][()]
    novel_replay_jumpiness = f['replay_jumpiness'][()]

novel_replay_ts = novel_replay_ts[novel_replay_jumpiness>0]

# Plot results
plt.figure(figsize=(3,1)) # TODO plot examples with Eric
plt.subplot(211)
plt.title('Pre-task REM (non-stationary events only)')
plt.imshow(
     REMpre_posteriors.T,
     aspect='auto',
     interpolation='none',
     cmap='magma',
     rasterized=True,
     vmin=0.02,
     vmax=.045
     )
plt.xlim(0,3600)
plt.xticks([0,1800,3600],['','',''])
plt.yticks([0,19,39],[0,50,100])
plt.ylabel('Location (cm)')

plt.subplot(413)
plt.eventplot(novel_replay_ts,colors='C2')
plt.yticks([])
plt.xlim(0,3600)
plt.xticks([0,1800,3600],[0,1,2])
plt.xlabel('Time (min)')
plt.ylabel('')

plt.savefig("../../output_REM/preplay_exampleBayesianReplay.pdf")

#%% Analysis
results_dir = params['path_to_output']+"/bayesian_replay"
resultsList = os.listdir(results_dir)

#%%
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
        replay_freqs = 1/(np.diff(h5_file['replay_locs'][()]/params['sampling_frequency']))
        replay_freqs = np.append(replay_freqs,0)
        for i in range(len(h5_file["replay_locs"][()])):
            data_list.append(  # This will create one list entry per cell
                {
                    "eventID": i,
                    "mouse": h5_file["mouse"][()].decode("utf-8"),
                    "condition": h5_file["condition"][()].decode("utf-8"),
                    "Type": "replay" if file_name.endswith("REMpost.h5") else "preplay",
                    "replayEventTime": h5_file['replay_locs'][i]/params['sampling_frequency'],
                    "replayFrequency": replay_freqs[i],
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

# %% Plot numEvents preplay vs replay
sns.boxenplot(
    data=df.query("replayEventJumpiness>0"),
    y='replayFrequency',
    x='Type',
    palette=['C3','C0'],
    showfliers=False
    #errorbar='se'
)
plt.title('Non-stationary only')
plt.xlabel('')
plt.yscale('log')
plt.ylabel('Replay frequency (Hz)')
plt.savefig("../../output_REM/REMpre_vs_post_replayFrequency.pdf")

#%% DESCRIPTIVES
mean_preplay = df.query("replayEventJumpiness>0 and Type=='preplay'")['replayFrequency'].mean()
SEM_preplay = df.query("replayEventJumpiness>0 and Type=='preplay'")['replayFrequency'].sem()
print(f'{mean_preplay} +/- {SEM_preplay}')

mean_replay = df.query("replayEventJumpiness>0 and Type=='replay'")['replayFrequency'].mean()
SEM_replay = df.query("replayEventJumpiness>0 and Type=='replay'")['replayFrequency'].sem()
print(f'{mean_replay} +/- {SEM_replay}')

#%% STATS
pg.ttest(
    x=df.query("replayEventJumpiness>0 and Type=='preplay'")['replayFrequency'],
    y=df.query("replayEventJumpiness>0 and Type=='replay'")['replayFrequency'],
    paired=True)

# %% Plot slope preplay vs replay
plt.figure(figsize=(.75,1))
sns.barplot(
    data=df_numEvents,
    y='numReplayEvents',
    x='Type',
    palette=['C3','C0'],
    #showfliers=False
    errorbar='se',
    capsize=.2
    
)
# plt.title('REM post')
plt.xlabel('')
plt.xticks([0,1],['Preplay', 'Replay'],rotation=90)
plt.ylabel('Num. flexible sequences')
plt.savefig("../../output_REM/REMpre_vs_post_numEvents.pdf")


#%%
sns.boxenplot(
    data=df.query("replayEventJumpiness>0"),
    y='replayEventSlope',
    x='Type',
    palette=['C3','C0'],
    showfliers=False
)
plt.title('REM post')
plt.ylabel('Replay event\nspeed (cm.s$^{-1}$)')
plt.xlabel('')
plt.savefig("../../output_REM/REMpre_vs_post_slope.pdf")

#%%
plt.figure(figsize=(1,.25))
sns.histplot(
    data=df.query("Type=='replay' and replayEventJumpiness>0"),
    x='replayEventScore',
)
# plt.title('REM post')
plt.xlabel('Replay score (R$^{2}$)')
plt.ylabel('N')
plt.savefig("../../output_REM/REMpost_replayScores.pdf")

#%% DESCRIPTIVES
mean_S1 = df.query("Type=='replay' and replayEventJumpiness>0")['replayEventScore'].mean()
SEM_S1 = df.query("Type=='replay' and replayEventJumpiness>0")['replayEventScore'].sem()
print(f'{mean_S1} +/- {SEM_S1}')

