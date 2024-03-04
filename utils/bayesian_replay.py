# %% Imports
import numpy as np
from numpy import polyfit

# %%
def extract_linear_replay(posterior_probs, params):
    replayLocs = []
    replayScore = []
    replayJumpiness = []
    replayLength = []

    # Get shuffled values
    for shuffle_i in range(params['numShuffles']):
        
        pass
        

    currentWindowIdx = np.array([0:params['windowSize']])
    while True:
        
        

        currentWindow = posterior_probs[currentWindowIdx]

        # For each window, compute score, jumpiness, length

        # If 3 scores exceed shuffled surrogate, append index and properties to variables
    
        current_window+=params['stepSize']
        if current_window[1]>len(posterior_probs):
            break

    return replayLocs, replayScore, replayJumpiness, replayLength

# %%