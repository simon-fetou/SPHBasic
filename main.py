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
from equations.timeStep import Euler_step
from Visualization.visu import init_plot, update_plot
from dataExtraction.extract import write_vtk

#assigning global variables for local use

lx_domain = prm.lx_domain
ly_domain = prm.ly_domain
Mx = prm.Mx
My = prm.My
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


def main():
    # Initialiation des tableaux
    pos = np.zeros((Npart, 2))
    vel = np.zeros((Npart, 2))
    accel = np.zeros((Npart, 2))
    drho_dt = np.zeros(Npart)
    drhou_dt = np.zeros((Npart,2))
    rho = np.zeros(Npart)
    pres = np.zeros(Npart)

    # Initialisation des particules
    pos, rho = init_part(Mx,My)
    # Question : est-ce que ipart est egal à Npart????

    #Saving the initial time vtk file
    write_vtk(f"champs{0}.vtk", pos, vel, accel, rho, pres)

    # Preparation de la figure pour l'animation
    plt.close("all")
    fig, cax1, cax2, cax3, scat1, scat2, scat3, quiv, time_text, time_template = init_plot()
    plt.show(block=False)
    plt.pause(0.1)

    # Boucle temporelle (avec affichage des resultats toutes les nsave iterations
    start = time.time()
    for n in range(1,N):
        pos,vel,rho,pres = Euler_step(n,pos,vel,rho,pres)

        #Saving the vtk file
        write_vtk(f"champs{n}.vtk", pos, vel, accel, rho, pres)

        if n%nsave==0:
            end = time.time()
            print('iteration {} time {:.5f} tpsCPU {:.2f}'.format(n,n*dt,(end-start)) )
            start = time.time()
            time_text.set_text(time_template%(n*dt))
            update_plot(pos, pres, rho, vel, fig, cax1, cax2, cax3, scat1, scat2, scat3, quiv)
            plt.draw()
            plt.pause(0.01)


if __name__ == "__main__":
    main()