import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from random import random as rand
import time

# Domaine et taille du patch initial
lx_patch =1.
ly_patch =1.
lx_domain = 2.
ly_domain = 2.

# Gravite
g = np.array([0., -0.5])

# Proprietes des fluides
rho_0 = 1. # densite de reference
pres0 = 1.
gamma = 1.
vel_ref = 0.5*np.sqrt(abs(g[1])*ly_patch)
c_ref = 10*vel_ref
k = rho_0*c_ref**2/gamma # coefficient dans l'equation de Tait

# Parametres numeriques
h = 0.1 # Smoothing length
alpha = 0.5  # Coefficient de viscosite artificielle
h_sur_dx = 0.9 # 

# Constantes du probleme
dx = h/0.9
dy = dx
Mx = int(lx_domain/dx)
My = int(ly_domain/dy)
Npart = int(Mx**2/2)
m_0 = rho_0*(dx**2)

# Parametres temporels
Tf = 4
dt = h/(8*c_ref)
N = int(Tf/dt)
nsave = 5