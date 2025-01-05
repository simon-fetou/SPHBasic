import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from random import random as rand
#from numba import jit
import time
import parameters as prm
from equations.kernel import length,gradW
from boundaryConditions.boundaryConditions import boundary

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


def Euler_step(n, pos, vel, rho, pres):

    drho_dt = np.zeros(Npart)
    drhou_dt = np.zeros((Npart,2))

    # Calcul des taux de variation de la densite drho_dt et de la QDM drhou_dt pour chaque particule i
    for i in range(Npart):
        drho_dt_i = float(0)
        drhou_dt_i = np.array([0., 0.])

        # on fait la sommation sur les voisins j
        for j in range(Npart):

            r_ij =pos[i]-pos[j] # vecteur position entre i et j            #r_ij=xi-xj
            d_ij = length(r_ij) # distance entre i et j

            if d_ij<=2*h: # test si la particule est dans le support compact 2h  (rij<=2h)

                assert (rho[j] != 0)  # test pour vérifier que rho_j n'est pas nul

                # Pressure : Tait equation
                pres[i] = k*((rho[i]/rho_0)**(gamma)-1)
                pres[j] = k*((rho[j]/rho_0)**(gamma)-1)

                # Equation d'Euler
                # taux de variation de la densite
                drho_dt_i += rho[i]*(vel[i]-vel[j]).dot(gradW(r_ij,h)*m_0/rho[j])

                # Calcul du terme de viscosite artificielle
                vij_rij =(vel[i]-vel[j]).dot(r_ij)/(d_ij**2) if (d_ij!=0) else 0 
                PIij = 0
                if vij_rij<0:
                    PIij =-alpha*h*c_ref*((rho[i]+rho[j])/2)*vij_rij
                # taux de variation de la QDM
                drhou_dt_i +=-(pres[i]+pres[j]+PIij)*(gradW(r_ij,h)*dx**2)#+rho[i]*g

        drho_dt[i] = drho_dt_i
        drhou_dt[i] = drhou_dt_i

    # Mise à jour de la densité, de la vitesse et de la position pour chaque particule i
    for i in range(Npart):
        rho[i] += drho_dt[i]*dt
        vel[i] += (drhou_dt[i]+g)*dt
        pos[i] += vel[i]*dt

        # Imposition des conditions aux limites (cf fonction plus bas)
        boundary(pos[i],vel[i])

    return pos, vel, rho, pres
