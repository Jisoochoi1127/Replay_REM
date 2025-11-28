# %% Imports
import numpy as np
from utils.bayesian_replay import linear_fit
from matplotlib import pyplot as plt

#%% Generate data
time_vec = np.arange(150)
position = np.random.randint(0,40,150)
position2 = np.arange(150)
position3 = np.arange(150,dtype='float')
idx_to_nan = np.random.choice(time_vec,20)
position3[idx_to_nan] = position3[idx_to_nan]*np.nan
# %%
plt.plot(time_vec,position3)
# %%
score, jumpiness, portion_window, slope = linear_fit(time_vec, position3)
# %%
posterior_probs=np.random.random((150,40))
MAP=np.argmax(posterior_probs,axis=1)
# %%
nan_col=np.ones(40)*np.nan
# %%
posterior_probs[6,:]=nan_col
posterior_probs[65,:]=nan_col
posterior_probs[42,:]=nan_col
# %%
nan_locs = np.isnan(np.max(posterior_probs,axis=1))
actual_map=np.max(posterior_probs,axis=1)
actual_map[nan_locs] = actual_map[nan_locs]*np.nan
# %%
score, jumpiness, portion_window, slope = linear_fit(time_vec, actual_map)
# %%
