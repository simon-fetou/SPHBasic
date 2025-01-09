import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import parameters as prm

#assigning global variables for local use

lxDomain = prm.lxDomain
lyDomain = prm.lyDomain
dx = prm.dx
dy = prm.dy

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


pos = np.zeros((Npart, 2))
vel = np.zeros((Npart, 2))
drhoDt = np.zeros(Npart)
drhoUDt = np.zeros((Npart,2))
rho = np.zeros(Npart)
press = np.zeros(Npart)


def init_plot():

    fig = plt.figure(figsize=(12,8))
    plt.set_cmap('jet')
    ax = fig.subplots(nrows=2, ncols=2, sharex=True, sharey=True)

    ax[0,0].set_xlim(0, lxDomain)
    ax[0,0].set_ylim(0, lyDomain)
    ax[0,0].set_xlabel('x')
    ax[0,0].set_ylabel('y')
    time_template = 'time = %.2fs'
    time_text = ax[0,0].text(0.5*lxDomain,1.1*lyDomain,'') 
    time_text.set_text(time_template%(0.))
    
    scat1 = ax[0,0].scatter(pos[:,0],pos[:,1],c=press)
    div = make_axes_locatable(ax[0,0])
    cax1 = div.append_axes('right', '5%', '5%')
    cax1 = fig.colorbar(scat1, ax=ax[0, 0], label='Pressure')

    scat2 = ax[0,1].scatter(pos[:,0],pos[:,1],c=rho)
    div = make_axes_locatable(ax[0,1])
    cax2 = div.append_axes('right', '5%', '5%')
    cax2 = fig.colorbar(scat2, ax=ax[0, 1], label='Density')

    quiv = ax[1,0].quiver(pos[:,0],pos[:,1], vel[:,0], vel[:,1],color=['r','b','g'], scale=21)

    scat3 = ax[1,1].scatter(pos[:,0],pos[:,1],c=np.sqrt(vel[:,0]**2+vel[:,1]**2))
    div = make_axes_locatable(ax[1,1])
    cax3 = div.append_axes('right', '5%', '5%')
    cax3 = fig.colorbar(scat3, ax=ax[1, 1], label='Velocity Magnitude')
    
 
    return fig, cax1, cax2, cax3, scat1, scat2, scat3, quiv, time_text, time_template


def update_plot(pos, pres, rho, vel, fig, cax1, cax2, cax3, scat1, scat2, scat3, quiv):

    scat1.set_offsets(pos)
    scat1.set_array(press)
    scat1.set_clim(vmin=np.min(press), vmax=np.max(press)) 
    #cax1.cla() # clear current axes (colormap)
    #fig.colorbar(scat1, cax=cax1,label='Pressure') # set new colormap with new values
    #scat1.autoscale() # Autoscale the scalar limits on the norm instance using the current array

    scat2.set_offsets(pos)
    scat2.set_array(rho)
    scat2.set_clim(vmin=np.min(rho), vmax=np.max(rho))
    #cax2.cla() # clear current axes (colormap)
    #fig.colorbar(scat2, cax=cax2,label='rho') # set new colormap with new values
    #scat2.autoscale() # Autoscale the scalar limits on the norm instance using the current array

    quiv.set_offsets(pos)
    quiv.set_UVC(vel[:,0],vel[:,1])
    #quiv.set_array(vel_mag)
    #quiv.autoscale() # Autoscalenorm = matplotlib.colors.Normalize()

    velMag = np.sqrt(vel[:,0]**2+vel[:,1]**2)

    scat3.set_offsets(pos)
    scat3.set_array(velMag)
    scat3.set_clim(vmin=np.min(velMag), vmax=np.max(velMag))
    #cax3.cla() # clear current axes (colormap)
    #fig.colorbar(scat3, cax=cax3,label='vel_mag') # set new colormap with new values
    #scat3.autoscale() # Autoscale the scalar limits on the norm instance using the current array

    return fig, cax1, cax2, cax3, scat1, scat2, scat3, quiv
