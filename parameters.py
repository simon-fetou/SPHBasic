import numpy as np

#-----------------------------Global parameters---------------------------
g = np.array([0., -9.81])                   # gravity

#------------------------Domain and Fluid patch parameters--------------------
x0D = 0.
y0D = 0.
lxDomain = 2.
lyDomain = 2.
lxFluid = lxDomain/2.
lyFluid = 1.75
x0F = 0.                                    #centering the fluid patch (x dir)
y0F = 0.                                    # from bottom of recepient (y dir)

#-----------------------------Kernel parameters---------------------------
h = 0.05                                     # Smoothing length
h_sur_dx = 1.4                              # number of fluid particles within h proximity of a particle

#-------------------------------Fluid properties----------------------------
rho0 = 1.                                   # reference density
press0 = 1.013                              # reference pressure
gamma = 1.
velRef = 0.5*np.sqrt(abs(g[1])*lyFluid)     # reference velocity
cRef = 10*velRef                            # reference sound velocity
k = rho0*cRef**2/gamma                      # TAIT equation coefficient
alpha = 0.5                                 # Artificial viscosity 

#-----------------------------Fluid patch discretization------------------
dx = h/h_sur_dx
dy = dx
Mx = int(lxFluid/dx)
My = int(lyFluid/dy)
Npart = int(Mx*My)
m0 = rho0*(dx*dy)                           #mass of a fluid particle

#-----------------------------Times parameters-----------------------------
Tf = 2
dt = h/(8*cRef)
Nt = int(Tf/dt)
nsave = 10                                  # saving interval