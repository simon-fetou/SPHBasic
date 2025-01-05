import numpy as np
import parameters as prm

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


def init_part(Mx,My):

    pos = np.zeros((Npart, 2))
    rho = np.zeros(Npart)
    ipart=0
    for i in range(int(Mx/2)):
        for j in range(My):
            x=(lx_domain/2-lx_patch/2)+i*dx+dx/2    #Pour centrer le patch 
            y=j*dx+dx/2
            pos[ipart]= x,y 
            rho[ipart]= rho_0*(1.+(gamma*g[1]/(c_ref**2)*(ly_patch-y)))**(1./gamma) #avec un abus sur la cond. hydrostatique
            ipart+=1
    return pos, rho