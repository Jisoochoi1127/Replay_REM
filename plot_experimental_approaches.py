# %%
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import yaml
import h5py
import scipy.io as sio
from utils.helperFunctions import load_data

plt.style.use("plot_style.mplstyle")

# %%
with open("params.yaml", "r") as file:
    params = yaml.full_load(file)

# %% Define conditions, dataset
states_list = ["REMpre", "wake", "REMpost"]
condition_list = ["LTD1", "LTD5", "HATD5"]
mouse_list = ["pv1060", "pv1254", "pv1069"]

# %% First, plot SFPs
mouse = "pv1060"
condition = "LTD1"
state = "wake"
data = load_data(mouse, condition, state, params)
# %%
plt.figure()
plt.imshow(np.max(data["SFPs"], axis=2), vmin=0, vmax=10, cmap="YlGnBu_r")
plt.axis("off")
plt.savefig("../../output_REM/SFPs.pdf")
# %% Then plot transients during each state
wake_data = load_data(mouse, condition, "wake", params)
REMpre_data = load_data(mouse, condition, "REMpre", params)
REMpost_data = load_data(mouse, condition, "REMpost", params)

# %% Extract ephys traces
with h5py.File('../../output_REM/LFP/lfp_pv1060_LTD1.h5','r') as f:
    pre_REM_LFP = f['pre_rem_lfp'][()]
    post_REM_LFP = f['post_rem_lfp'][()]

# %%
import matplotlib

# cmap = matplotlib.cm.get_cmap('nipy_spectral')
cmap = matplotlib.cm.get_cmap("viridis")
# numNeurons=data['SFPs'].shape[0]
numNeurons = 100
threshold = 0.2
LFP_Fs = 2000
MS_Fs = 30
interval_pre = [11,13]

plt.figure(figsize=(3, 1.5))
plt.subplot(434)
plt.title("REM pre", color='C1')
plt.plot(pre_REM_LFP[22000:26000,1],
         color = 'C1',
         linewidth=.3)
plt.axis("off")

plt.subplot(435)
plt.title("Wakefulness")
# animal position here?
plt.axis("off")

plt.subplot(436)
plt.title("REM post", color='C4')
plt.plot(post_REM_LFP[20000:24000,1],
         color = 'C4',
         linewidth=.3)
plt.axis("off")

plt.subplot(234)
plt.imshow(
    REMpre_data["binaryData"].T,
    cmap="bone_r",
    vmax=1.5,
    aspect="auto",
    interpolation="none"
)
plt.xlim(330, 390)
plt.xticks([330, 360, 390], [0, 1, 2])
plt.ylim(0, numNeurons)
# plt.xlabel('Time (s)')
plt.ylabel("Neuron #")


plt.subplot(235)
plt.imshow(wake_data["binaryData"].T,
           cmap="bone_r",
            vmax=1.5,
           aspect="auto",
           interpolation="none")

plt.xlim(0, 60)
plt.xticks([0, 30, 60], [0, 1, 2])
plt.ylim(0, numNeurons)
plt.yticks([])
plt.xlabel("Time (s)")
plt.ylabel("")


plt.subplot(236)
plt.imshow(
    REMpost_data["binaryData"].T,
    cmap="bone_r",
    vmax=1.5,
    aspect="auto",
    interpolation="none"
)

plt.xlim(300, 360)
plt.xticks([300, 330, 360], [0, 1, 2])
plt.ylim(0, numNeurons)
# plt.xlabel('Time (s)')
plt.yticks([])
plt.ylabel("")

plt.savefig("../../output_REM/calcium_transients.pdf")

# %% Plot anxiety track avoidance
HAT_data = {}
path = "../../datasets/REM_data/pv1060/HATD5"

f = sio.loadmat(path + "/behav.mat")
HAT_data.update(
    {  # Note that older files do not have background/tones
        "position": f["behav"]["position"][0][0],
        "behavTime": f["behav"]["time"][0][0].T[0] / 1000,
    }
)

HAT_data["LT_position"] = -1 * (HAT_data["position"][:, 0] - 100)

plt.figure(figsize=(1.5, 0.5))
plt.subplot(1, 2, 1)
plt.plot(HAT_data["behavTime"], HAT_data["LT_position"], color="C4")
plt.xticks([0, 300])
plt.xlabel("Time (s)")
plt.ylabel("Position\n(cm)")

plt.subplot(1, 4, 3)
sns.kdeplot(y=HAT_data["LT_position"], color="C4")
plt.xticks([0, 0.025], [0, 0.025])
plt.yticks([])
plt.savefig("../../output_REM/anxiety_behavior.pdf")

# %% Plot anxiety track avoidance
normal_data = {}
path = "../../datasets/REM_data/pv1060/LTD5"

f = sio.loadmat(path + "/behav.mat")
normal_data.update(
    {  # Note that older files do not have background/tones
        "position": f["behav"]["position"][0][0],
        "behavTime": f["behav"]["time"][0][0].T[0] / 1000,
    }
)

normal_data["LT_position"] = -1 * (normal_data["position"][:, 0] - 100)

plt.figure(figsize=(1.5, 0.5))
plt.subplot(1, 2, 1)
plt.plot(normal_data["behavTime"], normal_data["LT_position"])
plt.xticks([0, 300])
plt.xlabel("Time (s)")
plt.ylabel("Position\n(cm)")

plt.subplot(1, 4, 3)
sns.kdeplot(
    y=normal_data["LT_position"],
)
plt.xticks([0, 0.025], [0, 0.025])
plt.yticks([])
plt.savefig("../../output_REM/normal_behavior.pdf")
