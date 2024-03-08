# %% Imports
import numpy as np
from numpy import polyfit
#%%
def linear_fit(time_vec, position):
    # Assess replay length
    nanIdx = np.isnan(position)
    portion_window = sum(~nanIdx)/len(time_vec)

    # Ignore NaN values
    time_vec = time_vec[~nanIdx]
    position = position[~nanIdx]

    # Compute replay score
    _, residuals, _, _, _ = polyfit(x=time_vec, y=position, deg=1, full=True)
    score = abs(1 - residuals[0]/((len(position)-1) * np.var(position))) # Using absolute will consider both replay directions
    jumpiness = np.nanmax(np.diff(position)) # Maximum jump across two frames

    return score, jumpiness, portion_window

def extract_linear_replay(posterior_probs, params):
    np.random.seed(params['seed']) # For reproducibility

    replayLocs = []
    replayScore = []
    replayJumpiness = []
    replayPortion = []

    locationIdx = np.arange(posterior_probs.shape[1]) # Will be used to shuffle
    # Compute actual maximum a posteriori from posterior probabilities, convert to cm
    actual_map = (np.argmax(posterior_probs,axis=1)+params['spatialBinSize']/2)*params['spatialBinSize']

    # Shuffle posteriors, extract shuffled maps
    shuffled_maps = np.zeros((params['numShuffles'],len(posterior_probs)))*np.nan
    
    for shuffle_i in range(params['numShuffles']):
        # Shuffle posterior locations
        np.random.shuffle(locationIdx)
        shuffled_posteriors = posterior_probs[:,locationIdx]

        # Shuffle time
        split_point = np.random.randint(len(posterior_probs))
        shuffled_posteriors = np.concatenate((shuffled_posteriors[split_point:], shuffled_posteriors[:split_point]))
        
        # Compute argmax on shuffled data
        shuffled_maps[shuffle_i,:] = (np.argmax(shuffled_posteriors,axis=1)+params['spatialBinSize']/2)*params['spatialBinSize']

    ## Sliding window analysis
    # Initialize first window
    currentWindowIdx = np.arange(params['windowSize'])

    while currentWindowIdx[-1]<len(posterior_probs):
        # For each window, compute score, jumpiness, portion replayed
        actual_score, actual_jumpiness, actual_portion = linear_fit(currentWindowIdx, actual_map[currentWindowIdx])

        # Same for shuffled
        shuffled_score = np.zeros(params['numShuffles'])
        shuffled_jumpiness = np.zeros(params['numShuffles'])
        shuffled_portion = np.zeros(params['numShuffles'])
        for shuffle_i in range(params['numShuffles']):
            shuffled_score[shuffle_i], shuffled_jumpiness[shuffle_i], shuffled_portion[shuffle_i] = linear_fit(currentWindowIdx, shuffled_maps[shuffle_i,currentWindowIdx])
    
        # If 3 scores exceed shuffled surrogate, append index and properties to variables
        if actual_score>=np.percentile(shuffled_score, 95) and actual_jumpiness<=np.percentile(shuffled_jumpiness,5) and actual_portion>=np.percentile(shuffled_portion,95):
            if not replayLocs[-1]+params['windowSize'] > currentWindowIdx[0]: # This prevents counting the same event multiple times
                replayLocs.append(currentWindowIdx[0])
                replayScore.append(actual_score)
                replayJumpiness.append(actual_jumpiness)
                replayPortion.append(actual_portion)
        
        currentWindowIdx+=params['stepSize'] # Step forward

    return replayLocs, replayScore, replayJumpiness, replayPortion