from random import random as rand
import parameters as prm
from equations.kernel import length,gradW

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

def boundary(pos,vel):

    """
    Function sets the boundary conditions and particles behaviour (Dummy formulation to start)

    arguments:
    pos: [x,y] positions of fluid particles
    vel: [velx,vely] velocity of fluid particles

    Returns:
    pos:[x,y] updated positions of fluid particles 
    vel: [velx,vely] updated velocity of fluid particles
    """

    limit=0.5*h                  # closest limit to boundaries

    #Left wall
    if pos[0] - limit < 0:                  # if x+Little gets out of the domain
        
        if pos[0] - limit/2 < 0:                        #If x even closer to the Right bound
            pos[0] =0 + limit/2 + 0.01*rand()           # Send him back in

        vel[0] *= -0.5                      # Reorientate it's vel and dampen it by 50%

    #Right wall
    elif pos[0] + limit >= lxDomain:        # if x-Little gets out of the domain

        if pos[0] + limit/2 > lxDomain:                 #If x even closer to the Left bound
            pos[0] =lxDomain - (limit/2+0.01*rand())    # Send him back in

        vel[0] *= -0.5                      # Reorientate it's vel and dampen it by 50%
    
    #Bottom wall
    if pos[1] < limit:                      # if y+Little gets out of the domain

        if pos[1] -limit/2 < 0:                         #If y even closer to the Bottom bound
            pos[1] =0 + limit/2+0.01*rand()

        vel[1] *=-0.5                       # Reorientate it's vel and dampen it by 50%
    
    #Top wall 
    elif pos[1] + limit >= lyDomain:        # if y-Little gets out of the domain

        if pos[1] + limit/2 > lyDomain:                 #If y even closer to the Top bound
            pos[1] =lyDomain - (limit/2+0.01*rand())

        vel[1] *= -0.5                      # Reorientate it's vel and dampen it by 50%

    return pos,vel