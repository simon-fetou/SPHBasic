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

lxDomain = prm.lxDomain
lyDomain = prm.lyDomain
dx = prm.dx
dy = prm.dy

x0F = prm.x0F
y0F = prm.y0F
lxFluid = prm.lxFluid
lyFluid = prm.lyFluid

h = prm.h
alpha = prm.alpha

cRef = prm.cRef

m0 = prm.m0
k = prm.k

rho0 = prm.rho0

dt = prm.dt
Tf = prm.Tf
Nt = prm.Nt
nsave= prm.nsave

gamma = prm.gamma
g = prm.g
Npart = prm.Npart


def EulerTimeStep(n, pos, vel, rho, press):
    """
    Function solves the evolution of quantities from time n to n+1 with an Euler approach

    Arguments:
    pos: [[x,y],...] positions of fluid particles (n)
    vel: [[Vx,Vy],...] velocities of fluid particles (n)
    rho: [rhoi,...] densities of fluid particles (n)
    press: [pressi,...] pressures of fluid particles (n)

    Returns:
    pos: [[x,y],...] positions of fluid particles (n+1)
    vel: [[Vx,Vy],...] velocities of fluid particles (n+1)
    rho: [rhoi,...] densities of fluid particles (n+1)
    press: [pressi,...] pressures of fluid particles (n+1)
    """

    drhoDt = np.zeros(Npart)
    drhoUDt = np.zeros((Npart,2))

    # Density and momentum variation over time (drhoDt,drhoUDt) for each particle i (NS equations)
    for i in range(Npart):
        drho_dt_i = float(0)
        drhou_dt_i = np.array([0., 0.])

        # Summing up over all neighbor particles j (i included)
        for j in range(Npart):

            r_ij =pos[i]-pos[j]                 # position vector between i and j        
            d_ij = length(r_ij)                 # distance between i et j

            if d_ij<=2*h: # ensuring j is within the compact bounds 2h (rij<=2h kernel formulation)
                          # Helps limit computing cost

                assert (rho[j] != 0)            # rho_j should never be null (divisions)

                # Pressure : Tait equation
                press[i] = k*((rho[i]/rho0)**(gamma)-1)
                press[j] = k*((rho[j]/rho0)**(gamma)-1)

                # Euler Equation
                    # density variation
                drho_dt_i += rho[i]*(vel[i]-vel[j]).dot(gradW(r_ij,h)*m0/rho[j])

                    # Artificial viscosity
                vij_rij =(vel[i]-vel[j]).dot(r_ij)/(d_ij**2) if (d_ij!=0) else 0 
                PIij = 0

                if vij_rij<0:
                    PIij =-alpha*h*cRef*((rho[i]+rho[j])/2)*vij_rij

                    # Momentum variation
                drhou_dt_i +=-(press[i]+press[j]+PIij)*(gradW(r_ij,h)*dx**2)#+rho[i]*g

        drhoDt[i] = drho_dt_i
        drhoUDt[i] = drhou_dt_i

    # Updating quantities
    for i in range(Npart):
        rho[i] += drhoDt[i]*dt
        vel[i] += (drhoUDt[i]+g)*dt
        pos[i] += vel[i]*dt

        # Enforcing boundary conditions
        boundary(pos[i],vel[i])

    return pos, vel, rho, press
