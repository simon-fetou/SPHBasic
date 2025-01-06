import numpy as np
import parameters as prm
from otherFunctions.NewtonRaphson import newtonRaphson_poly

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
press0 = prm.press0

dt = prm.dt
Tf = prm.Tf
Nt = prm.Nt
nsave= prm.nsave

gamma = prm.gamma
g = prm.g
Npart = prm.Npart


def init_part(Mx,My):
    """
    function initializes fluid particles in the domain 

    Arguments:
    Mx: Number of fluid particles in x direction
    My: Number of fluid particles in y direction
    
    Returns:
    pos: pos[i] = [x,y] position of the fluid particles
    rho: rho[i] = [rho] density of the fluid particles
    press: press[i] = [press] pressure of the fluid particles
    """

    pos = np.zeros((Npart, 2))
    rho = np.zeros(Npart)
    press = np.zeros(Npart)
    ipart=0

    for i in range(int(Mx)):
        for j in range(My):

            x= x0F + i*dx + dx/2        # x-coordinate of Fluid patch from left start 
            y= y0F + j*dy + dy/2        # y-coordinate of Fluid patch from bottom start
            pos[ipart]= x,y 

            #initializing density of particle respecting hydrostatic and TAIT equations
                #coefficients (ax^n + bx + c = 0)
            a = 1
            b = -(gamma*(rho0**(gamma-1))*np.abs(g[1])/(cRef**2))*(lyFluid-pos[ipart][1])
            c= -(rho0**gamma)*(1 + gamma*press0/(rho0*(cRef**2)))
            n= gamma
                #Root guess
            rhoG = rho0                                          

            rho[ipart]=newtonRaphson_poly(a, b, c, n, rhoG)

            press[ipart] = press0 + rho[ipart]*np.abs(g[1])*(lyFluid-pos[ipart][1])

            ipart+=1

    return pos, rho, press