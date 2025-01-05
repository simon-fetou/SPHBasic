from random import random as rand
import parameters as prm
from equations.kernel import length,gradW

#assigning global variables for local use

lx_domain = prm.lx_domain
ly_domain = prm.ly_domain
dx = prm.dx
dy = prm.dy

lx_patch = prm.lx_patch
ly_patch = prm.ly_patch

h = prm.h
alpha = prm.alpha

c_ref = prm.c_ref

m_0 = prm.m_0
k = prm.k

rho_0 = prm.rho_0

dt = prm.dt
Tf = prm.Tf
N = prm.N
nsave= prm.nsave

gamma = prm.gamma
g = prm.g
Npart = prm.Npart

def boundary(pos,vel):
    limit=0.5*h
    if pos[0] >=lx_domain-limit:  # si x >= fin domaine - limite
        
        if pos[0] <lx_domain-limit/2: # si x < fin domaine - limite/2
            pos[0] =lx_domain - (limit/2+0.01*rand()) # position x bloquee a fin domaine  - (limite/2 + 0.01*rand())
        vel[0] *= -0.5   # la composante x de la vitesse est reflechie et amortie de 50%
    elif pos[0] <limit: # debut domaine
        if pos[0] <limit/2:
            pos[0] =0 + limit/2+0.01*rand()
        vel[0] *= -0.5

    if pos[1] >=ly_domain-limit:
        if pos[1] <ly_domain-limit/2:
            pos[1] =ly_domain - (limit/2+0.01*rand())
        vel[1] *=-0.5
    elif pos[1] < limit:
        if pos[1] < limit/2:
            pos[1] =0 + limit/2+0.01*rand()
        vel[1] *= -0.5
    return pos,vel