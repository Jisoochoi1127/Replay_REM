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

def extract_linear_replay_shuffle_types(posterior_probs, params):
    np.random.seed(params['seed']) # For reproducibility
    replayLocs_P = []
    replayScore_P = []
    replayJumpiness_P = []
    replayPortion_P = []
    replaySlope_P = []
    replayLocs_T = []
    replayScore_T = []
    replayJumpiness_T = []
    replayPortion_T = []
    replaySlope_T = []
    replayLocs_PT = [] 
    replayScore_PT = []
    replayJumpiness_PT = []
    replayPortion_PT = []
    replaySlope_PT = []

    positionIdx = np.arange(posterior_probs.shape[0])
    locationIdx = np.arange(posterior_probs.shape[1]) # Will be used to shuffle
    # Compute actual maximum a posteriori from posterior probabilities, convert to cm
    actual_map = (np.argmax(posterior_probs,axis=1)+params['spatialBinSize']/2)*params['spatialBinSize']

    # Shuffle posteriors, extract shuffled maps
    shuffled_maps_T = np.zeros((params['numShuffles'],len(posterior_probs)))*np.nan #shuffle time
    shuffled_maps_P = np.zeros((params['numShuffles'],len(posterior_probs)))*np.nan #shuffle position
    shuffled_maps_PT = np.zeros((params['numShuffles'],len(posterior_probs)))*np.nan #shuffle time and position
    
    for shuffle_i in range(params['numShuffles']):
        # Shuffle position
        np.random.shuffle(locationIdx)
        np.random.shuffle(positionIdx)
        shuffled_posteriors_T = posterior_probs[:,locationIdx] # time only
        shuffled_posteriors_P = posterior_probs[positionIdx,:] # position only
        shuffled_posteriors_PT = posterior_probs[positionIdx,:] # position...
        shuffled_posteriors_PT = shuffled_posteriors_PT[:,locationIdx] # ... and time
        
        # Compute argmax on shuffled data
        shuffled_maps_P[shuffle_i,:] = (np.argmax(shuffled_posteriors_P,axis=1)+params['spatialBinSize']/2)*params['spatialBinSize']
        shuffled_maps_T[shuffle_i,:] = (np.argmax(shuffled_posteriors_T,axis=1)+params['spatialBinSize']/2)*params['spatialBinSize']
        shuffled_maps_PT[shuffle_i,:] = (np.argmax(shuffled_posteriors_PT,axis=1)+params['spatialBinSize']/2)*params['spatialBinSize']

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
        shuffled_score_P = np.zeros(params['numShuffles']) 
        shuffled_jumpiness_P = np.zeros(params['numShuffles']) 
        shuffled_portion_P = np.zeros(params['numShuffles']) 
        shuffled_score_T = np.zeros(params['numShuffles']) 
        shuffled_jumpiness_T = np.zeros(params['numShuffles']) 
        shuffled_portion_T = np.zeros(params['numShuffles']) 
        shuffled_score_PT = np.zeros(params['numShuffles']) 
        shuffled_jumpiness_PT = np.zeros(params['numShuffles']) 
        shuffled_portion_PT = np.zeros(params['numShuffles'])

        for shuffle_i in range(params['numShuffles']):
            shuffled_score_P[shuffle_i], shuffled_jumpiness_P[shuffle_i], shuffled_portion_P[shuffle_i], _ = linear_fit(
                currentWindowIdx/params['sampling_frequency'],
                shuffled_maps_P[shuffle_i,currentWindowIdx]
                )
            shuffled_score_T[shuffle_i], shuffled_jumpiness_T[shuffle_i], shuffled_portion_T[shuffle_i], _ = linear_fit(
                currentWindowIdx/params['sampling_frequency'],
                shuffled_maps_T[shuffle_i,currentWindowIdx]
                )
            shuffled_score_PT[shuffle_i], shuffled_jumpiness_PT[shuffle_i], shuffled_portion_PT[shuffle_i], _ = linear_fit(
                currentWindowIdx/params['sampling_frequency'],
                shuffled_maps_PT[shuffle_i,currentWindowIdx]
                )
    
        # If scores and jumpiness exceed shuffled surrogate, append index and properties to variables
        # FOR POSITION
        if actual_score>=np.percentile(shuffled_score_P, 95) and actual_jumpiness<=np.percentile(shuffled_jumpiness_P,5) and actual_portion>=np.percentile(shuffled_portion_P,95):
            if not replayLocs_P or replayLocs_P[-1]+params['windowSize'] <= currentWindowIdx[0]:
                replayLocs_P.append(currentWindowIdx[0])
                replayScore_P.append(actual_score)
                replayJumpiness_P.append(actual_jumpiness)
                replayPortion_P.append(actual_portion)
                replaySlope_P.append(actual_slope)

        # FOR TIME
        if actual_score>=np.percentile(shuffled_score_T, 95) and actual_jumpiness<=np.percentile(shuffled_jumpiness_T,5) and actual_portion>=np.percentile(shuffled_portion_T,95):
            if not replayLocs_T or replayLocs_T[-1]+params['windowSize'] <= currentWindowIdx[0]:
                replayLocs_T.append(currentWindowIdx[0])
                replayScore_T.append(actual_score)
                replayJumpiness_T.append(actual_jumpiness)
                replayPortion_T.append(actual_portion)
                replaySlope_T.append(actual_slope)

        # FOR POSITION AND TIME
        if actual_score>=np.percentile(shuffled_score_PT, 95) and actual_jumpiness<=np.percentile(shuffled_jumpiness_PT,5) and actual_portion>=np.percentile(shuffled_portion_PT,95):
            if not replayLocs_PT or replayLocs_PT[-1]+params['windowSize'] <= currentWindowIdx[0]:
                replayLocs_PT.append(currentWindowIdx[0])
                replayScore_PT.append(actual_score)
                replayJumpiness_PT.append(actual_jumpiness)
                replayPortion_PT.append(actual_portion)
                replaySlope_PT.append(actual_slope)
        
        currentWindowIdx+=params['stepSize'] # Step forward

    output_dict = {
        'position_shuffle': {
            'replayLocs_P': replayLocs_P, 
            'replayScore_P':replayScore_P,
            'replayJumpiness_P': replayJumpiness_P,
            'replayPortion_P': replayPortion_P,
            'replaySlope_P': replaySlope_P
        },
        'time_shuffle': {
            'replayLocs_T': replayLocs_T, 
            'replayScore_T': replayScore_T,
            'replayJumpiness_T': replayJumpiness_T,
            'replayPortion_T': replayPortion_T,
            'replaySlope_T': replaySlope_T
        },
        'position_time_shuffle': {
            'replayLocs_PT': replayLocs_PT, 
            'replayScore_PT': replayScore_PT,
            'replayJumpiness_PT': replayJumpiness_PT,
            'replayPortion_PT': replayPortion_PT,
            'replaySlope_PT': replaySlope_PT
        }
    }

    return output_dict

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