import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from random import random as rand
import time
import parameters as prm
from boundaryConditions.boundaryConditions import boundary
from initConditions.initConditions import init_part
from equations.kernel import length, gradW
from equations.timeStep import EulerTimeStep
from Visualization.visu import init_plot, update_plot
from dataExtraction.extract import write_vtk

#assigning global variables for local use

lxDomain = prm.lxDomain
lyDomain = prm.lyDomain
Mx = prm.Mx
My = prm.My
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


def main():
    # Initialiation des tableaux
    pos = np.zeros((Npart, 2))
    vel = np.zeros((Npart, 2))
    accel = np.zeros((Npart, 2))
    drhoDt = np.zeros(Npart)
    drhoUDt = np.zeros((Npart,2))
    rho = np.zeros(Npart)
    press = np.zeros(Npart)

    # Initialisation des particules
    pos, rho, pres = init_part(Mx,My)
    # Question : est-ce que ipart est egal à Npart????

    #Saving the initial time vtk file
    write_vtk(f"champs{0}.vtk", pos, vel, accel, rho, press)

    # Boucle temporelle (avec affichage des resultats toutes les nsave iterations
    start = time.time()
    for n in range(1,Nt):

        pos,vel,rho,pres = EulerTimeStep(n,pos,vel,rho,press)

        if n%nsave==0:

            #Saving the vtk file
            write_vtk(f"champs{n}.vtk", pos, vel, accel, rho, press)

            end = time.time()
            print('iteration {} time {:.5f} tpsCPU {:.2f}'.format(n,n*dt,(end-start)) )
            start = time.time()

if __name__ == "__main__":
    main()