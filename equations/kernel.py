import numpy as np

def length(r):
    return np.sqrt(r.dot(r))

def gradW(r, h):              # gradient Cubic kernel 2H
    C=30/(14*np.pi*h**2)      #on est en 2D
    q=length(r)/h
    gradW=C*r/(h**2)*((-2+3*q/2)*(q<1)+((-(0.5/(q+0.001*(q==0)))*(2-q)**2))*((q>1)*(q<2)))
    return gradW
