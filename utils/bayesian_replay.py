# %% Imports
import numpy as np
from numpy import polyfit, polyval

# %% Polyfit analysis
def linear_fit(time_vec, position):

    # Assess replay length
    nanIdx = np.isnan(position)
    portion_window = sum(~nanIdx)/len(time_vec)

    # Ignore NaN values
    time_vec = time_vec[~nanIdx]
    position = position[~nanIdx]

    # Compute replay score
    _, residuals, _, _, _ = polyfit(x=time_vec, y=position, deg=1, full=True)
    score = 1 - residuals[0]/((len(position)-1) * np.var(position))
    jumpiness = np.nanmax(np.diff(position)) # Maximum jump across two frames

    return score, jumpiness, portion_window

# %%
def extract_linear_replay(posterior_probs, params):
    replayLocs = []
    replayScore = []
    replayJumpiness = []
    replayLength = []

    #TODO deal with preallocation, empty windows

    currentWindowIdx = np.arange(params['windowSize'])

    while True:
        currentWindow = posterior_probs[currentWindowIdx]

        # Compute maximum a posteriori, output results in cm for correct speed computations
        #TODO assert empty/nan frames remain nan values
        actual_map = (np.argmax(currentWindow,axis=1)+params['spatialBinSize']/2)*params['spatialBinSize']

        # For each window, compute score, jumpiness, length
        actual_score, actual_jumpiness, actual_length = linear_fit(currentWindowIdx, actual_map)

        
    
        currentWindowIdx+=params['stepSize'] # Step forward

        if currentWindowIdx[1]>len(posterior_probs): # If window goes beyond recording
            break

    # Get shuffled values
    for shuffle_i in range(params['numShuffles']):
        #TODO shuffle data here
        
        currentWindowIdx = np.array([0:params['windowSize']])
        while True:
            currentWindow = posterior_probs[currentWindowIdx]

            actual_map = (np.argmax(currentWindow,axis=1)+params['spatialBinSize']/2)*params['spatialBinSize']

            # For each window, compute score, jumpiness, length
            shuffled_score, shuffled_jumpiness, shuffled_length = linear_fit(time,actual_map)


            current_window+=params['stepSize']
            if current_window[1]>len(posterior_probs):
                break
        
    # If 3 scores exceed shuffled surrogate, append index and properties to variables
    significant_replay_events = np.zeros(len(actual_score))
    for i, time_loc in enumerate(significant_replay_events):
        if actual_score>=np.percentile(shuffled_score,95) and actual_jumpiness<np.percentile(shuffled_jumpiness,95) and actual_length>np.percentile(shuffled_length,95)
            replayLocs.append(#TODO)
            replayScore.append(actual_score[i])
            replayJumpiness.append(actual_jumpiness[i])
    return replayLocs, replayScore, replayJumpiness, replayLength

