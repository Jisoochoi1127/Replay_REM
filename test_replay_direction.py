#%%
import numpy as np
from numpy import polyfit
import matplotlib.pyplot as plt

#%%
time_vec=np.arange(15)
position=np.arange(15)*3
reversed_position = np.flip(position)

#%%
_, residuals, _, _, _ = polyfit(x=time_vec, y=position, deg=1, full=True)
normal_score = 1 - residuals[0]/((len(position)-1) * np.var(position))
reversed_score = 1 - residuals[0]/((len(position)-1) * np.var(position))
print('R$^{2}_{normal}$ = '+ f'{normal_score}')
print('R$^{2}_{reversed}$ = '+ f'{reversed_score}')
# %%
