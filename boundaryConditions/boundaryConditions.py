import numpy as np
import math
from random import random as rand
import parameters as prm
from equations.kernel import length,gradW

#assigning global variables for local use

x0D = prm.x0D
y0D = prm.y0D
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
press0 = prm.press0

dt = prm.dt
Tf = prm.Tf
Nt = prm.Nt
nsave= prm.nsave

gamma = prm.gamma
g = prm.g
Npart = prm.Npart


def boundary(pos,vel,rho,press):

    '''
    Function sets the boundary conditions and particles behaviour (Dummy formulation to start)

    arguments:
    pos: [x,y] positions of fluid particles
    vel: [velx,vely] velocity of fluid particles

    Returns:
    pos:[x,y] updated positions of fluid particles 
    vel: [velx,vely] updated velocity of fluid particles
    '''

    limit=0.5*h                  # closest limit to boundaries

    #Left wall
    if pos[0] - limit < x0D:                     # if x+Little gets out of the domain
        
        leftWall = [x0D , pos[1]]              # takes the vector between particle and wall
        doi = np.linalg.norm(pos-leftWall)     # distance from the wall
        toi = doi/h

        rho = rho * (1-0.5*math.erfc(toi))
        press = press0 + k*((rho/rho0)**(gamma)-1)
        vel[0] *= -0.5                      # Dampen it's vel more and more till 0 on the B (non ph)

    #Right wall
    elif pos[0] + limit >= lxDomain:           # if x-Little gets out of the domain

        rightWall = [lxDomain , pos[1]]      # takes the vector between particle and wall
        doi = np.linalg.norm(pos-rightWall)
        toi = doi/h

        rho = rho * (1-0.5*math.erfc(toi))
        press = press0 + k*((rho/rho0)**(gamma)-1)
        vel[0] *= -0.5                      # Dampen it's vel more and more till 0 on the B (non ph)
    
    #Bottom wall
    if pos[1] - limit < y0D:                      # if y+Little gets out of the domain

        bottomWall = [pos[0] , y0D]      # takes the vector between particle and wall
        doi = np.linalg.norm(pos-bottomWall)
        toi = doi/h

        rho = rho * (1-0.5*math.erfc(toi))
        press = press0 + k*((rho/rho0)**(gamma)-1)
        vel[1] *= -0.5                    # Dampen it's vel more and more till 0 on the B (non ph)
    
    #Top wall 
    elif pos[1] + limit >= lyDomain:           # if y-Little gets out of the domain

        topWall = [pos[0] , lyDomain]          # takes the vector between particle and wall
        doi = np.linalg.norm(pos-topWall)      # distance from the wall
        toi = doi/h

        rho = rho * (1-0.5*math.erfc(toi))
        press = press0 + k*((rho/rho0)**(gamma)-1)
        vel[1] *= -0.5                    # Dampen it's vel more and more till 0 on the B (non ph)

    return pos,vel,rho,press