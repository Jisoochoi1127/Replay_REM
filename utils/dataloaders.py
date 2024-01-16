import h5py
import numpy as np
import scipy.io as sio
import os

def load_data(path):
    data = {}
    
    # Basic information
    data.update(
                {
                'path': path,
                'task': 'LT'
                }
    )

    try: # If recent MATLAB format
        f = h5py.File(path + '/ms.mat','r')
        data.update({
                    'SFPs':np.array(f.get('ms/SFPs')),
                    'caTime':np.array(f.get('ms/time'))[0]/1000, # convert ms->s
                    'rawData':np.array(f.get('ms/RawTraces')).T})
    except: # If legacy MATLAB format
        f = sio.loadmat(path + '/ms.mat')
        data.update(
                    {
                    'SFPs':f['ms']['SFPs'][0][0],
                    'caTime':f['ms']['time'][0][0].T[0]/1000, # convert ms->s
                    'rawData':f['ms']['RawTraces'][0][0] 
                    }
                    )

    try: # If recent MATLAB format
        f = h5py.File(path + '/behav.mat','r')
        data.update(
                    {
                    'position':np.array(f.get('behav/position')).T,
                    'behavTime':np.array(f.get('behav/time'))[0]/1000, # convert ms->s
                    }
                    )

    except:
        f = sio.loadmat(path + '/behav.mat')
        try: # If old format
            data.update(
                        { # Note that older files do not have background/tones
                        'position':f['behav']['position'][0][0],
                        'behavTime':f['behav']['time'][0][0].T[0]/1000,
                        }
                        )
        except: # else must be recent Deeplabcut output
            data.update(
                        { # Note that older files do not have background/tones
                        'position':f['behav']['ledPosition'][0][0], # Use LED to match older recordings
                        }
                        )
            if len(f['behav']['time'][0][0][0])>1:
                data.update(
                        {
                            'behavTime':f['behav']['time'][0][0][0]/1000
                         }
                         )
            else:
                data.update(
                        {
                            'behavTime':f['behav']['time'][0][0].T[0]/1000
                        }
                )

    return data