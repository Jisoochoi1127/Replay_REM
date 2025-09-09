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
    score = np.nan
    jumpiness = np.nan
    slope = np.nan
    if len(position)>0:
        p, residuals, _, _, _ = polyfit(x=time_vec, y=position, deg=1, full=True)
        if len(residuals)>0:
            score = 1 - residuals[0]/((len(position)-1) * np.var(position))
            jumpiness = np.nanmax(np.diff(position)) # Maximum jump across two frames
            slope = p[0]
        
    return score, jumpiness, portion_window, slope

def extract_linear_replay(posterior_probs, params):
    np.random.seed(params['seed']) # For reproducibility

    replayLocs = []
    replayScore = []
    replayJumpiness = []
    replayPortion = []
    replaySlope = []

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
        actual_score, actual_jumpiness, actual_portion, actual_slope = linear_fit(
            currentWindowIdx/params['sampling_frequency'], # Divide to compute slope in cm/s
            actual_map[currentWindowIdx]
            )

        # Same for shuffled
        shuffled_score = np.zeros(params['numShuffles'])
        shuffled_jumpiness = np.zeros(params['numShuffles'])
        shuffled_portion = np.zeros(params['numShuffles'])
        for shuffle_i in range(params['numShuffles']):
            shuffled_score[shuffle_i], shuffled_jumpiness[shuffle_i], shuffled_portion[shuffle_i], _ = linear_fit(
                currentWindowIdx/params['sampling_frequency'],
                shuffled_maps[shuffle_i,currentWindowIdx]
                )
    
        # If scores and jumpiness exceed shuffled surrogate, append index and properties to variables
        if actual_score>=np.percentile(shuffled_score, 95) and actual_jumpiness<=np.percentile(shuffled_jumpiness,5) and actual_portion>=np.percentile(shuffled_portion,95):
            if not replayLocs or replayLocs[-1]+params['windowSize'] <= currentWindowIdx[0]:
                replayLocs.append(currentWindowIdx[0])
                replayScore.append(actual_score)
                replayJumpiness.append(actual_jumpiness)
                replayPortion.append(actual_portion)
                replaySlope.append(actual_slope)
        
        currentWindowIdx+=params['stepSize'] # Step forward

    return replayLocs, replayScore, replayJumpiness, replayPortion, replaySlope


def extract_raw_linear_replay(sorted_binary, params):
    np.random.seed(params['seed']) # For reproducibility

    replayLocs = []
    replayScore = []
    replayJumpiness = []
    replayPortion = []
    replaySlope = []
    
    # Find the identity of neuron currently being active
    active_neuron = np.argmax(sorted_binary,axis=1).astype('float') #TODO if multiple neurons active, use their average identity?

    # Remove silent epochs
    active_neuron[np.sum(sorted_binary,axis=1)==0]=np.nan

    shuffled_active_neuron = np.zeros((params['numShuffles'],len(sorted_binary)))*np.nan

    for shuffle_i in range(params['numShuffles']):
        # Shuffle cell identities as recommended in Foster, 2017
        cellIdx = np.arange(sorted_binary.shape[1])
        np.random.shuffle(cellIdx)
        shuffled_binary = sorted_binary[:,cellIdx]
        
        # Compute argmax on shuffled data
        shuffled_active_neuron[shuffle_i,:] = np.argmax(shuffled_binary,axis=1).astype('float') # idetify active neuron
        shuffled_active_neuron[shuffle_i,np.sum(shuffled_binary,axis=1)==0] = np.nan # remove silent epochs

    ## Sliding window analysis
    # Initialize first window
    currentWindowIdx = np.arange(params['windowSize'])

    while currentWindowIdx[-1]<len(sorted_binary):
        shuffled_maps = np.zeros((params['numShuffles'],len(sorted_binary)))*np.nan
        # For each window, compute score, jumpiness, portion replayed

        actual_score, actual_jumpiness, actual_portion, actual_slope = linear_fit(
            currentWindowIdx/params['sampling_frequency'], # Divide to compute slope in cm/s
            active_neuron[currentWindowIdx]
            )

        # Same for shuffled
        shuffled_score = np.zeros(params['numShuffles'])
        shuffled_jumpiness = np.zeros(params['numShuffles'])
        shuffled_portion = np.zeros(params['numShuffles'])
        for shuffle_i in range(params['numShuffles']):
            shuffled_score[shuffle_i], shuffled_jumpiness[shuffle_i], shuffled_portion[shuffle_i], _ = linear_fit(
                currentWindowIdx/params['sampling_frequency'],
                shuffled_active_neuron[shuffle_i,currentWindowIdx]
                )
    
        # If scores and jumpiness exceed shuffled surrogate, append index and properties to variables
        if actual_score>=np.percentile(shuffled_score, 95) and actual_jumpiness<=np.percentile(shuffled_jumpiness,5) and actual_portion>=np.percentile(shuffled_portion,95):
            if not replayLocs or replayLocs[-1]+params['windowSize'] <= currentWindowIdx[0]:
                replayLocs.append(currentWindowIdx[0])
                replayScore.append(actual_score)
                replayJumpiness.append(actual_jumpiness)
                replayPortion.append(actual_portion)
                replaySlope.append(actual_slope)
        
        currentWindowIdx+=params['stepSize'] # Step forward

    return replayLocs, replayScore, replayJumpiness, replayPortion, replaySlope