import numpy as np
import math
from random import random as rand
import parameters as prm
from equations.kernel import length,W, gradW

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


def reflect_particle(pos:np.ndarray, wall_position:float, axis:float):
    '''
    Function reflects the position of a particle close to the wall (to build ghost particles).

    Parameters:
    pos: fluid particle position [x, y].
    wall_position: Position of the wall along the given axis.
    axis: 0 for x-axis, 1 for y-axis.

    Return: 
    Reflected position of the ghost particle.
    '''

    reflected_position = np.copy(pos)
    reflected_position[axis] = 2 * wall_position - pos[axis]         # demonstrable

    return reflected_position


def boundary(pos:np.ndarray, vel:np.ndarray, rho:np.ndarray, press:np.ndarray):

    '''
    Function sets the boundary conditions and particles behaviour

    arguments:
    pos: [x,y] positions of fluid particles
    vel: [velx,vely] velocity of fluid particles
    dens: [dens] density of fluid particles
    press: [press] pressure of fluid particles

    Returns:
    Np_ghost: number of ghost particles
    ghost_pos: [x,y] positions of ghost particles
    vel: [velx,vely] velocity of ghost particles
    dens: [dens] density of ghost particles
    press: [press] pressure of ghost particles

    '''
    
    ghost_pos = np.empty((0,2), dtype=float)
    ghost_vel = np.empty((0,2), dtype=float)
    ghost_rho =np.array([])
    ghost_press = np.array([])

    limit = 0.5*h
    #Left wall
    if pos[0] - limit < x0D:                                # if x-Little gets out of the domain

        #Updating density, and pressure of the particle at the boundary
        leftWall = [x0D , pos[1]]                           # takes the LeftWall position
        doi = np.linalg.norm(pos-leftWall)                  # distance from the wall
        toi = doi/h
        rho = rho * (1-0.5*math.erfc(toi))
        press = press0 + k*((rho/rho0)**(gamma)-1)

        # creating a ghost particle to update particles velocities on the boundaries later
        reflected_pos = reflect_particle(pos, x0D, 0)
        ghost_pos =np.vstack([ghost_pos, reflected_pos])
        ghost_vel = np.vstack([ghost_vel, -vel])                            #no-slip wall
        ghost_rho = np.append(ghost_rho, [rho])
        ghost_press = np.append(ghost_press, [press])
        
    #Right wall
    elif pos[0] + limit >= lxDomain:           # if x+Little gets out of the domain

        #Updating density, and pressure of the particle at the boundary
        rightWall = [lxDomain , pos[1]]                           # takes the RightWall position
        doi = np.linalg.norm(pos-rightWall)
        toi = doi/h
        rho = rho * (1-0.5*math.erfc(toi))
        press = press0 + k*((rho/rho0)**(gamma)-1)

        # creating ghost a particle to update particles velocities on the boundaries later
        reflected_pos = reflect_particle(pos, lxDomain, 0)
        ghost_pos =np.vstack([ghost_pos, reflected_pos])
        ghost_vel = np.vstack([ghost_vel, -vel])                            #no-slip wall
        ghost_rho = np.append(ghost_rho, [rho])
        ghost_press = np.append(ghost_press, [press])

    #Bottom wall
    if pos[1] - limit < y0D:                      # if y-Little gets out of the domain

        #Updating density, and pressure of the particle at the boundary
        bottomWall = [pos[0] , y0D]                           # takes the BottomWall position
        doi = np.linalg.norm(pos-bottomWall)
        toi = doi/h
        rho = rho * (1-0.5*math.erfc(toi))
        press = press0 + k*((rho/rho0)**(gamma)-1)

        # creating a ghost particle to update particles velocities on the boundaries later
        reflected_pos = reflect_particle(pos, y0D, 1)
        ghost_pos =np.vstack([ghost_pos, reflected_pos])
        ghost_vel = np.vstack([ghost_vel, -vel])                            #no-slip wall
        ghost_rho = np.append(ghost_rho, [rho])
        ghost_press = np.append(ghost_press, [press])
        
    #Top wall 
    elif pos[1] + limit >= lyDomain:           # if y+Little gets out of the domain
        
        #Updating density, and pressure of the particle at the boundary
        topWall = [pos[0] , lyDomain]                           # takes the RightWall position
        doi = np.linalg.norm(pos-topWall)                       # distance from the wall
        toi = doi/h
        rho = rho * (1-0.5*math.erfc(toi))
        press = press0 + k*((rho/rho0)**(gamma)-1)

        # creating a ghost particle to update particles velocities on the boundaries later
        reflected_pos = reflect_particle(pos, lyDomain, 1)
        ghost_pos =np.vstack([ghost_pos, reflected_pos])
        ghost_vel = np.vstack([ghost_vel, -vel])                            #no-slip wall
        ghost_rho = np.append(ghost_rho, [rho])
        ghost_press = np.append(ghost_press, [press])

    # Get the number of ghost particles Np_ghost
    Np_ghost = len(ghost_pos)

    velGh = 0.0
    for j in range(Np_ghost):
        r_ij = pos-ghost_pos[j]            # position vector between i and the ghost j 
        d_ij = length(r_ij)                # distance between i and ghost j

        if d_ij<=2*h:     # ensuring ghost j is within the compact bounds 2h (rij<=2h kernel formulation)
                          # Helps limit computing cost

            assert (ghost_rho[j] != 0)            # rho_j should never be null (divisions)

            velGh += ghost_vel[j]*W(r_ij,h)*m0/ghost_rho[j]


    vel += velGh


    return pos,vel,rho,press
