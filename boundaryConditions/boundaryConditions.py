import numpy as np
import math
from random import random as rand
import parameters as prm
from equations.kernel import length,W, gradW
from dataExtraction.extract import write_vtk

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


def reflectAboutaxis(pos:np.ndarray, wall_position:float, axis:float):
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

def reflectAboutPoint(pos:np.ndarray, point:np.ndarray):
    '''
    Function reflects the position of a particle in walls angles (to build ghost particles).
    The particle is close to two walls

    Parameters:
    pos: fluid particle position [x, y].
    point: Position of the angle point.
    
    Return: 
    Reflected position of the ghost particle.
    '''

    reflected_position = np.copy(pos)
    reflected_position[0] = pos[0] +2*(point[0]-pos[0])         # demonstrable
    reflected_position[1] = pos[1] +2*(point[1]-pos[1])         # demonstrable

    return reflected_position


def ghostMaker(pos:np.ndarray, vel:np.ndarray, rho:np.ndarray, press:np.ndarray):

    '''
    Function builds the  ghost particles in solid boundaries

    arguments:
    pos: [[x,y],...] positions of fluid particles
    vel: [[velx,vely],...] velocity of fluid particles
    rho: [[rho],...] density of fluid particles
    press: [[press],...] pressure of fluid particles

    Returns:
    Np_ghost: number of ghost particles
    ghost_pos: [[x,y],...] positions of ghost particles
    ghost_vel: [[velx,vely],...] velocity of ghost particles
    ghost_rho: [[rho],...] density of ghost particles
    ghost_press: [[press],...] pressure of ghost particles
    
    '''

    ghost_pos = np.empty((0,2), dtype=float)
    ghost_vel = np.empty((0,2), dtype=float)
    ghost_rho =np.array([])
    ghost_press = np.array([])

    limit = 2*h        # making sure even a particle on a wall can build a compact kernel with 2h 
                       #smooth length 

    for i in range(Npart):
        
        #-----------------------------Wall about ghost particles----------------------------------
        #Left wall
        if (pos[i][0] - limit < x0D):                        # if x-Little gets out of the domain
        
            # creating a ghost particle to update particles velocities close to boundaries later
            reflected_pos = reflectAboutaxis(pos[i], x0D, 0)
            ghost_pos =np.vstack([ghost_pos, reflected_pos])
            ghost_vel = np.vstack([ghost_vel, [-vel[i][0],-vel[i][1]]])         #no-slip wall condition
            ghost_rho = np.append(ghost_rho, [rho[i]])
            ghost_press = np.append(ghost_press, [press[i]])
            
        #Right wall
        if (pos[i][0] + limit > lxDomain):           # if x+Little gets out of the domain

            # creating ghost a particle to update particles velocities on the boundaries later
            reflected_pos = reflectAboutaxis(pos[i], lxDomain, 0)
            ghost_pos =np.vstack([ghost_pos, reflected_pos])
            ghost_vel = np.vstack([ghost_vel, [-vel[i][0],-vel[i][1]]])            #no-slip wall
            ghost_rho = np.append(ghost_rho, [rho[i]])
            ghost_press = np.append(ghost_press, [press[i]])
            

        #Bottom wall
        if (pos[i][1] - limit < y0D):                      # if y-Little gets out of the domain
        
            # creating a ghost particle to update particles velocities on the boundaries later
            reflected_pos = reflectAboutaxis(pos[i], y0D, 1)
            ghost_pos =np.vstack([ghost_pos, reflected_pos])
            ghost_vel = np.vstack([ghost_vel, [-vel[i][0],-vel[i][1]]])             #no-slip wall
            ghost_rho = np.append(ghost_rho, [rho[i]])
            ghost_press = np.append(ghost_press, [press[i]])
            
        #Top wall 
        if (pos[i][1] + limit > lyDomain):           # if y+Little gets out of the domain
        
            # creating a ghost particle to update particles velocities on the boundaries later
            reflected_pos = reflectAboutaxis(pos[i], lyDomain, 1)
            ghost_pos = np.vstack([ghost_pos, reflected_pos])
            ghost_vel = np.vstack([ghost_vel, [-vel[i][0],-vel[i][1]]])             #no-slip wall
            ghost_rho = np.append(ghost_rho, [rho[i]])
            ghost_press = np.append(ghost_press, [press[i]])
            
        #-----------------------------Angle about ghost particles----------------------------------
        #Bottom left angle point
        if (pos[i][0] - limit < x0D) and (pos[i][1] - limit < y0D):
            anglePoint = [x0D,y0D]
            reflected_pos = reflectAboutPoint(pos[i],anglePoint)
            ghost_pos = np.vstack([ghost_pos, reflected_pos])
            ghost_vel = np.vstack([ghost_vel, [-vel[i][0],-vel[i][1]]])             #no-slip wall
            ghost_rho = np.append(ghost_rho, [rho[i]])
            ghost_press = np.append(ghost_press, [press[i]])
            
        #Bottom right angle point
        if (pos[i][0] + limit > lxDomain) and (pos[i][1] - limit < y0D):
            anglePoint = [lxDomain,y0D]
            reflected_pos = reflectAboutPoint(pos[i],anglePoint)
            ghost_pos = np.vstack([ghost_pos, reflected_pos])
            ghost_vel = np.vstack([ghost_vel, [-vel[i][0],-vel[i][1]]])             #no-slip wall
            ghost_rho = np.append(ghost_rho, [rho[i]])
            ghost_press = np.append(ghost_press, [press[i]])
            
        #Top left angle point
        if (pos[i][0] - limit < x0D) and (pos[i][1] + limit > lyDomain):
            anglePoint = [x0D,lyDomain]
            reflected_pos = reflectAboutPoint(pos[i],anglePoint)
            ghost_pos = np.vstack([ghost_pos, reflected_pos])
            ghost_vel = np.vstack([ghost_vel, [-vel[i][0],-vel[i][1]]])              #no-slip wall
            ghost_rho = np.append(ghost_rho, [rho[i]])
            ghost_press = np.append(ghost_press, [press[i]])
            
        #Top right angle point
        if (pos[i][0] + limit > lxDomain) and (pos[i][1] + limit > lyDomain):
            anglePoint = [lxDomain,lyDomain]
            reflected_pos = reflectAboutPoint(pos[i],anglePoint)
            ghost_pos = np.vstack([ghost_pos, reflected_pos])
            ghost_vel = np.vstack([ghost_vel, [-vel[i][0],-vel[i][1]]])               #no-slip wall
            ghost_rho = np.append(ghost_rho, [rho[i]])
            ghost_press = np.append(ghost_press, [press[i]])
            
    # Get the number of ghost particles Np_ghost
    Np_ghost = len(ghost_pos)

    return Np_ghost,ghost_pos,ghost_vel,ghost_rho,ghost_press

"""
def boundary(pos:np.ndarray, vel:np.ndarray, rho:np.ndarray, press:np.ndarray):

    '''
    Function sets the boundary conditions and particles behaviour

    arguments:
    pos: [x,y] positions of fluid particles
    vel: [velx,vely] velocity of fluid particles
    dens: [dens] density of fluid particles
    press: [press] pressure of fluid particles

    Returns:
    pos: [x,y] positions of fluid particles
    vel: [velx,vely] velocity of fluid particles
    dens: [dens] density of fluid particles
    press: [press] pressure of fluid particles

    '''
    
    ghost_pos = np.empty((0,2), dtype=float)
    ghost_vel = np.empty((0,2), dtype=float)
    ghost_rho =np.array([])
    ghost_press = np.array([])
    
    limit = 0.5*h
    #Left wall
    if pos[0] - limit < x0D:                                # if x-Little gets out of the domain
        
        # creating a ghost particle to update particles velocities on the boundaries later
        reflected_pos = reflect_particle(pos, x0D, 0)
        ghost_pos =np.vstack([ghost_pos, reflected_pos])
        ghost_vel = np.vstack([ghost_vel, -vel])                            #no-slip wall
        ghost_rho = np.append(ghost_rho, [rho])
        ghost_press = np.append(ghost_press, [press])
        
    #Right wall
    elif pos[0] + limit >= lxDomain:           # if x+Little gets out of the domain

        # creating ghost a particle to update particles velocities on the boundaries later
        reflected_pos = reflect_particle(pos, lxDomain, 0)
        ghost_pos =np.vstack([ghost_pos, reflected_pos])
        ghost_vel = np.vstack([ghost_vel, -vel])                            #no-slip wall
        ghost_rho = np.append(ghost_rho, [rho])
        ghost_press = np.append(ghost_press, [press])

    #Bottom wall
    if pos[1] - limit < y0D:                      # if y-Little gets out of the domain
        
        # creating a ghost particle to update particles velocities on the boundaries later
        reflected_pos = reflect_particle(pos, y0D, 1)
        ghost_pos =np.vstack([ghost_pos, reflected_pos])
        ghost_vel = np.vstack([ghost_vel, -vel])                            #no-slip wall
        ghost_rho = np.append(ghost_rho, [rho])
        ghost_press = np.append(ghost_press, [press])
        
    #Top wall 
    elif pos[1] + limit >= lyDomain:           # if y+Little gets out of the domain
        
        # creating a ghost particle to update particles velocities on the boundaries later
        reflected_pos = reflect_particle(pos, lyDomain, 1)
        ghost_pos =np.vstack([ghost_pos, reflected_pos])
        ghost_vel = np.vstack([ghost_vel, -vel])                            #no-slip wall
        ghost_rho = np.append(ghost_rho, [rho])
        ghost_press = np.append(ghost_press, [press])

    # Get the number of ghost particles Np_ghost
    Np_ghost = len(ghost_pos)

    print('Np_ghost', Np_ghost)

    velGh = np.array([0.,0.])
    rhoGh = float(0.)
    pressGh = float(0.)

    for j in range(Np_ghost):
        r_ij = pos-ghost_pos[j]            # position vector between i and the ghost j 
        d_ij = length(r_ij)                # distance between i and ghost j

        if d_ij<=2*h:     # ensuring ghost j is within the compact bounds 2h (rij<=2h kernel formulation)
                          # Helps limit computing cost

            assert (ghost_rho[j] != 0)            # rho_j should never be null (divisions)

            velGh += ghost_vel[j]*W(r_ij,h)*m0/ghost_rho[j]
            rhoGh += ghost_rho[j]*W(r_ij,h)*m0/ghost_rho[j]
            pressGh += ghost_press[j]*W(r_ij,h)*m0/ghost_rho[j]


    vel += velGh
    rho += rhoGh
    press += pressGh


    return pos,vel,rho,press
"""

def boundary(pos:np.ndarray, vel:np.ndarray, rho:np.ndarray, press:np.ndarray):

    '''
    Function sets the boundary conditions and particles behaviour at the boundaries

    arguments:
    pos: [[x,y],...] positions of fluid particles
    vel: [[velx,vely],...] velocity of fluid particles
    rho: [[rho],...] density of fluid particles
    press: [[press],...] pressure of fluid particles

    Returns:
    pos: [[x,y],...] positions of fluid particles
    vel: [[velx,vely],...] velocity of fluid particles
    rho: [[rho],...] density of fluid particles
    press: [[press],...] pressure of fluid particles

    '''

    Np_ghost,ghost_pos,ghost_vel,ghost_rho,ghost_press = ghostMaker(pos, vel, rho, press)
    
    limit = 2*h        # From here close to boundaries the ghost particles should be taken into account

    for i in range(Npart):

        velGh = np.array([0.,0.])
        rhoGh = float(0.)
        pressGh = float(0.)

        if (pos[i][0] - limit < x0D) or (pos[i][0] + limit > lxDomain)\
           or (pos[i][1] - limit < y0D) or (pos[i][1] + limit > lyDomain):
            
            for j in range(Np_ghost):
                r_ij = pos[i]-ghost_pos[j]            # position vector between i and the ghost j 
                d_ij = length(r_ij)                # distance between i and ghost j

                if d_ij<=2*h:     # ensuring ghost j is within the compact bounds 2h (rij<=2h kernel formulation)
                          # Helps limit computing cost

                    assert (ghost_rho[j] != 0)            # rho_j should never be null (divisions)

                    velGh += ghost_vel[j]*W(r_ij,h)*m0/ghost_rho[j]
                    rhoGh += ghost_rho[j]*W(r_ij,h)*m0/ghost_rho[j]
                    pressGh += ghost_press[j]*W(r_ij,h)*m0/ghost_rho[j]


            vel[i] += velGh
            rho[i] += rhoGh
            press[i] += pressGh
    
    '''
    Np_ghost,ghost_pos,ghost_vel,ghost_rho,ghost_press = ghostMaker(pos, vel, rho, press)
    limit = 0.5*h

    for i in range(Npart):
        #Left wall
        if pos[i][0] - limit < x0D:                                # if x-Little gets out of the domain

            if pos[i][0] - limit/2 < x0D:
                pos[i][0] = pos[i][0] + limit/2 + 0.01*rand()
        
            vel[i][0] *=-0.5
        
        #Right wall
        elif pos[i][0] + limit >= lxDomain:           # if x+Little gets out of the domain

            if pos[i][0] + limit/2 >= lxDomain:
                pos[i][0] = pos[i][0] - limit - 0.01*rand()
        
            vel[i][0] *=-0.5

        #Bottom wall
        if pos[i][1] - limit < y0D:                      # if y-Little gets out of the domain
        
            if pos[i][1] - limit/2 < y0D:
                pos[i][1] = pos[i][1] = limit + 0.01*rand()
        
            vel[i][1] *=-0.5
        
        #Top wall 
        elif pos[i][1] + limit >= lyDomain:           # if y+Little gets out of the domain

            if pos[i][1] + limit/2 >= lyDomain:
                pos[i][1] = pos[i][1] - limit - 0.01*rand()
            
            vel[i][1] *=-0.5
        '''
    
    return pos,ghost_pos,vel,rho,press


"""
def boundary(Np_ghost:float, pos:np.ndarray, vel:np.ndarray, rho:np.ndarray, press:np.ndarray,
            ghost_pos:np.ndarray, ghost_vel:np.ndarray, ghost_rho:np.ndarray, ghost_press:np.ndarray):

    '''
    Function sets the boundary conditions and particles behaviour

    arguments:
    pos: [x,y] positions of fluid particles
    vel: [velx,vely] velocity of fluid particles
    dens: [dens] density of fluid particles
    press: [press] pressure of fluid particles

    Returns:
    pos: [x,y] positions of fluid particles
    vel: [velx,vely] velocity of fluid particles
    dens: [dens] density of fluid particles
    press: [press] pressure of fluid particles

    '''
    
    
    limit = 0.5*h

    velGh = np.array([0.,0.])
    rhoGh = float(0.)
    pressGh = float(0.)

    for j in range(Np_ghost):
        r_ij = pos-ghost_pos[j]            # position vector between i and the ghost j 
        d_ij = length(r_ij)                # distance between i and ghost j

        if d_ij<=2*h:     # ensuring ghost j is within the compact bounds 2h (rij<=2h kernel formulation)
                          # Helps limit computing cost

            assert (ghost_rho[j] != 0)            # rho_j should never be null (divisions)

            velGh += ghost_vel[j]*W(r_ij,h)*m0/ghost_rho[j]
            rhoGh += ghost_rho[j]*W(r_ij,h)*m0/ghost_rho[j]
            pressGh += ghost_press[j]*W(r_ij,h)*m0/ghost_rho[j]


    vel += velGh
    rho += rhoGh
    press += pressGh


    return pos,vel,rho,press

"""