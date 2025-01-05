import numpy as np
import parameters as prm
from otherFunctions.NewtonRaphson import newtonRaphson_poly

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
pres0 = prm.pres0

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
    pres = np.zeros(Npart)
    ipart=0

    for i in range(int(Mx/2)):
        for j in range(My):
            
            x=(lx_domain/2-lx_patch/2)+i*dx+dx/2    #Pour centrer le patch 
            y=j*dx+dx/2
            pos[ipart]= x,y 

            #initializing density of all particles respecting hydrostatic and TAIT equations
            a = 1                      #coefficient
            b = -(gamma*(rho_0**(gamma-1))*np.abs(g[1])/(c_ref**2))*(ly_patch-pos[ipart][1])
            c= -(rho_0**gamma)*(1 + gamma*pres0/(rho_0*(c_ref**2)))
            n= gamma
            rhoG = rho_0                                           # root guess
            rho[ipart]=newtonRaphson_poly(a, b, c, n, rhoG)

            pres[ipart] = pres0 + rho[ipart]*np.abs(g[1])*(ly_patch-pos[i][1])

            #rho[ipart]=newtonRaphson_poly()
            #rho[ipart]= rho_0*(1.+(gamma*g[1]/(c_ref**2)*(ly_patch-y)))**(1./gamma) #avec un abus sur la cond. hydrostatique
            ipart+=1

    return pos, rho, pres