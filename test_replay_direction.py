#%%
import numpy as np
from numpy import polyfit
import matplotlib.pyplot as plt

#%%
time_vec=np.arange(15)
position=np.arange(15)*3
flat=np.ones(15)+20
reversed_position = np.flip(position)

#%%
p, residuals, _, _, _ = polyfit(x=time_vec, y=position, deg=1, full=True)
normal_score = 1 - residuals[0]/((len(position)-1) * np.var(position))

p_r, residuals_r, _, _, _ = polyfit(x=time_vec, y=reversed_position, deg=1, full=True)
reversed_score = 1 - residuals_r[0]/((len(reversed_position)-1) * np.var(reversed_position))

p_f, residuals_f, _, _, _ = polyfit(x=time_vec, y=flat, deg=1, full=True)
flat_score = 1 - residuals_f[0]/((len(flat)-1) * np.var(flat))

print('R$^{2}_{normal}$ = '+ f'{normal_score}')
print('R$^{2}_{reversed}$ = '+ f'{reversed_score}')
print('R$^{2}_{flat}$ = '+ f'{flat_score}')
# %%
plt.plot(time_vec, position, label='right')
plt.plot(time_vec, reversed_position, label='left')
plt.plot(time_vec, flat, label='left')
# %%
