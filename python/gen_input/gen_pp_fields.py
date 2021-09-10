import numpy as np
import h5py

#---------------------------------------------------
datadir = "/pauli-storage/tblommel/He-VB2PP/paper_recreate/tau0/"
filename = "PPinp_nt4000dt01_break.h5"
#---------------------------------------------------

#---------------------------------------------------
time_scale = np.longfloat(2.4188843265857e-17)
eV_scale = np.longfloat( 27.211386245988)
Efield_scale = np.longfloat(5.14220674763e11)
#---------------------------------------------------

#---------------------------------------------------
Delta_P = np.longfloat(15e-15) / time_scale
omega_P = np.longfloat(0.57) / eV_scale
E_0     = np.longfloat(6.6e9) / Efield_scale
Delta_p = np.longfloat(0.5e-15) / time_scale
omega_p = np.longfloat(22) / eV_scale
e_0     = np.longfloat(8.6e7) / Efield_scale
dt      = 0.01
tau     = np.longfloat(0e-15) / time_scale
#---------------------------------------------------

nt = 4000
delay = 5
Efield = np.zeros((3, nt+1+delay))
efield = np.zeros((3, nt+1+delay))

f = h5py.File(datadir + filename, 'w')

for T in np.arange(nt+1+delay):
  t = dt * (T-delay)
  if(0 <= t <= Delta_P):
    Efield[0,T] = E_0 * np.sin(np.pi * t / Delta_P)**2 * np.sin(omega_P * t)
  if(tau <= t <= tau + Delta_p): 
    efield[0,T] = e_0 * np.sin(np.pi * (t-tau) / Delta_p)**2 * np.sin(omega_p * (t-tau))
  
f['E'] = Efield
f['e'] = efield
