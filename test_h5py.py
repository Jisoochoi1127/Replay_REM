#%%
import h5py
import numpy as np

#%%
data_dict = {
    0:[0,65,2,6,77,4234],
    1:[3,5,1523,1,523,5364]
}
# %%
with h5py.File('../../test.h5','w') as f:
    dict_group = f.create_group('timestamps')
    for k, v in data_dict.items():
        dict_group[str(k)] = v
    
# %%
f = h5py.File('../../test.h5','f')


# %%
