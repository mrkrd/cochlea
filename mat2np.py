import numpy as np
import scipy.io

mat_lin = scipy.io.loadmat('param_lin.mat')
mat_res = scipy.io.loadmat('param_res.mat')

np.savez('bm_pars.npz',
         Ls=mat_lin['Ls'].astype(float),
         Rs=mat_lin['Rs'].astype(float),
         Ct=mat_lin['Ct'].astype(float),
         Rbm=mat_lin['Rbm'].astype(float),
         Cbm=mat_lin['Cbm'].astype(float),
         Lbm=mat_lin['Lbm'].astype(float),
         Rh=mat_lin['Rh'].astype(float),
         Lh=mat_lin['Lh'].astype(float),
         ampl_corr=mat_lin['ampl_corr'].astype(float),
         Abm=mat_lin['Abm'].astype(float),
         Cbm=mat_lin['Cbm'].astype(float),
         freq_map=mat_res['freq_map_res'].astype(float),
         Qmin=mat_res['Qmin'].astype(float),
         Qmax=mat_res['Qmax'].astype(float),
         SAT1=mat_res['SAT1'].astype(float),
         SAT4=mat_res['SAT4'].astype(float))


