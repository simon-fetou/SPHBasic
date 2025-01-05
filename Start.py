import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from random import random as rand
#from numba import jit
import time 

#%matplotlib tk
"""
*****_______****
*****_______****
"""

def init_plot():
    fig = plt.figure(figsize=(12,8))
    plt.set_cmap('jet')
    ax = fig.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
    ax[0,0].set_xlim(0, lx_domain)
    ax[0,0].set_ylim(0, ly_domain)
    ax[0,0].set_xlabel('x')
    ax[0,0].set_ylabel('y')
    time_template = 'time = %.2fs'
    time_text = ax[0,0].text(0.5*lx_domain,1.1*ly_domain,'') 
    time_text.set_text(time_template%(0.))
    scat1 = ax[0,0].scatter(pos[:,0],pos[:,1],c=pres)
    div = make_axes_locatable(ax[0,0])
    cax1 = div.append_axes('right', '5%', '5%')

    scat2 = ax[0,1].scatter(pos[:,0],pos[:,1],c=rho)
    div = make_axes_locatable(ax[0,1])
    cax2 = div.append_axes('right', '5%', '5%')

    quiv = ax[1,0].quiver(pos[:,0],pos[:,1], vel[:,0], vel[:,1], color=['r','b','g'], scale=21)

    scat3 = ax[1,1].scatter(pos[:,0],pos[:,1],c=np.sqrt(vel[:,0]**2+vel[:,1]**2))
    div = make_axes_locatable(ax[1,1])
    cax3 = div.append_axes('right', '5%', '5%')
 
    return fig, cax1, cax2, cax3, scat1, scat2, scat3, quiv, time_text, time_template


def update_plot(pos, pres, rho, vel):
    scat1.set_offsets(pos)
    scat1.set_array(pres)
    cax1.cla() # clear current axes (colormap)
    fig.colorbar(scat1, cax=cax1,label='Pressure') # set new colormap with new values
    #scat1.autoscale() # Autoscale the scalar limits on the norm instance using the current array

    scat2.set_offsets(pos)
    scat2.set_array(rho)
    cax2.cla() # clear current axes (colormap)
    fig.colorbar(scat2, cax=cax2,label='rho') # set new colormap with new values
    #scat2.autoscale() # Autoscale the scalar limits on the norm instance using the current array

    quiv.set_offsets(pos)
    quiv.set_UVC(vel[:,0],vel[:,1])
    vel_mag = np.sqrt(vel[:,0]**2+vel[:,1]**2)
    quiv.set_array(vel_mag)
    #quiv.autoscale() # Autoscalenorm = matplotlib.colors.Normalize()

    scat3.set_offsets(pos)
    scat3.set_array(vel_mag)
    cax3.cla() # clear current axes (colormap)
    fig.colorbar(scat3, cax=cax3,label='vel_mag') # set new colormap with new values
    #scat3.autoscale() # Autoscale the scalar limits on the norm instance using the current array

    return fig, cax1, cax2, cax3, scat1, scat2, scat3, quiv

"""
*****_______****
*****_______****
"""

def length(r):
    return np.sqrt(r.dot(r))

def gradW(r, h):              # gradient Cubic kernel 2H
    C=30/(14*np.pi*h**2)      #on est en 2D
    q=length(r)/h
    gradW=C*r/(h**2)*((-2+3*q/2)*(q<1)+((-(0.5/(q+0.001*(q==0)))*(2-q)**2))*((q>1)*(q<2)))
    return gradW

"""
*****_______****
*****_______****
"""

def Euler_step(n, pos, vel, rho, pres):

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
                vij_rij =(vel[i]-vel[j]).dot(r_ij)/(d_ij**2)*(d_ij!=0)
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

"""
*****_______****
*****_______****
"""
def boundary(pos,vel):
    limit=0.5*h
    if pos[0] >=lx_domain-limit:  # si x >= fin domaine - limite
        
        if pos[0] <lx_domain-limit/2: # si x < fin domaine - limite/2
            pos[0] =lx_domain - (limit/2+0.01*rand()) # position x bloquee a fin domaine  - (limite/2 + 0.01*rand())
        vel[0] *= -0.5   # la composante x de la vitesse est reflechie et amortie de 50%
    elif pos[0] <limit: # debut domaine
        if pos[0] <limit/2:
            pos[0] =0 + limit/2+0.01*rand()
        vel[0] *= -0.5

    if pos[1] >=ly_domain-limit:
        if pos[1] <ly_domain-limit/2:
            pos[1] =ly_domain - (limit/2+0.01*rand())
        vel[1] *=-0.5
    elif pos[1] < limit:
        if pos[1] < limit/2:
            pos[1] =0 + limit/2+0.01*rand()
        vel[1] *= -0.5
    return pos,vel

"""
*****_______****
*****_______****
"""

def init_part(Mx,My):
    ipart=0
    for i in range(int(Mx/2)):
        for j in range(My):
            x=(lx_domain/2-lx_patch/2)+i*dx+dx/2    #Pour centrer le patch 
            y=j*dx+dx/2
            pos[ipart]= x,y 
            rho[ipart]= rho_0*(1.+(gamma*g[1]/(c_ref**2)*(ly_patch-y)))**(1./gamma) #avec un abus sur la cond. hydrostatique
            ipart+=1
    return pos, rho

"""
*****_______****
*****_______****
"""

#Domaine et taille du patch initial
lx_patch =1.
ly_patch =1.
lx_domain = 2.
ly_domain = 2.

# Gravite
g = np.array([0., -0.5])

# Proprietes des fluides
rho_0 = 1. # densite de reference
gamma = 1.
vel_ref = 0.5*np.sqrt(abs(g[1])*ly_patch)
c_ref = 10*vel_ref
k = rho_0*c_ref**2/gamma # coefficient dans l'equation de Tait

# Parametres numeriques
h = 0.1 # Smoothing length
alpha = 0.5  # Coefficient de viscosite artificielle
h_sur_dx = 0.9 # 

# Constantes du probleme
dx = h/0.9
Mx = int(lx_domain/dx)
My = int(ly_domain/dx)
Npart = int(Mx**2/2)
m_0 = rho_0*(dx**2)

# Parametres temporels
Tf = 2
dt = h/(8*c_ref)
N = int(Tf/dt)
nsave = 1

# Initialiation des tableaux
pos = np.zeros((Npart, 2))
vel = np.zeros((Npart, 2))
drho_dt = np.zeros(Npart)
drhou_dt = np.zeros((Npart,2))
rho = np.zeros(Npart)
pres = np.zeros(Npart)

# Initialisation des particules
pos, rho = init_part(Mx,My)
# Question : est-ce que ipart est egal à Npart????


# Preparation de la figure pour l'animation
plt.close("all")
fig, cax1, cax2, cax3, scat1, scat2, scat3, quiv, time_text, time_template = init_plot()
plt.show(block=False)
plt.pause(0.1)

# Boucle temporelle (avec affichage des resultats toutes les nsave iterations
start = time.time()
for n in range(1,N):
    pos,vel,rho,pres = Euler_step(n,pos,vel,rho,pres)
    #print('iteration',n,'time',n*dt)
    if n%nsave==0:
        end = time.time()
        print('iteration {} time {:.5f} tpsCPU {:.2f}'.format(n,n*dt,(end-start)) )
        start = time.time()
        time_text.set_text(time_template%(n*dt))
        update_plot(pos, pres, rho, vel)
        plt.draw()
        plt.pause(0.01)