# %% Imports
import numpy as np
from numpy import polyfit

# %%
def extract_linear_replay(posterior_probs, params):
    replayLocs = []
    replayScore = []
    replayJumpiness = []
    replayLength = []
    time = np.arange(params['windowSize'])

    #TODO deal with preallocation, empty windows

    currentWindowIdx = np.array([0:params['windowSize']])
    while True:
        currentWindow = posterior_probs[currentWindowIdx]

        # Compute maximum a posteriori, output results in cm for correct speed computations
        #TODO ensure empty/nan frames remain nan values
        actual_map = (np.argmax(currentWindow,axis=1)+params['spatialBinSize']/2)*params['spatialBinSize']

        # For each window, compute score, jumpiness, length
        actual_score, actual_jumpiness, actual_length = linear_fit(time,actual_map)

        
    
        current_window+=params['stepSize']
        if current_window[1]>len(posterior_probs):
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

# %% Polyfit analysis
def linear_fit(time, position):

# x=time;
# y=decoding.REM_decoded_position(time);

# %Ignoring nan values
# nanIndex=isnan(y);
# y(nanIndex) = [];
# x = x(~nanIndex);

# % check the smapling point
# sampling_percentage(window)=length(y)/length(time);


# %calculate Rsquare
# p = polyfit(x,y,1);

# slope(window)=p(1);
# intercept(window)=p(2);
# yfit = polyval(p,x);
# yresid = y - yfit;
# SSresid = sum(yresid.^2);
# SStotal = (length(y)-1) * var(y);
# rsq(window) = 1 - SSresid/SStotal;
# Replay_score_actual=rsq;\
    return score, jumpiness, length