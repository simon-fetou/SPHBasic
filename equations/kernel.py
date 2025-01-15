import numpy as np

def length(r):
    return np.sqrt(r.dot(r))

def W(r,h):
    
    #facteur de normalisation
    C=30/(14*np.pi*(h**2))    
    #paramètre q
    q=length(r)/h
    
    W=C*((2./3-q**2+(q**3)/2)*(q>0)*(q<1)+(((2-q)**3)/6)*((q>1)*(q<2))+0*(q>2))
    
    return W

def gradW(r, h):              # gradient Cubic kernel 2H
    C=30/(14*np.pi*h**2)      #on est en 2D
    q=length(r)/h
    gradW=C*r/(h**2)*((-2+3*q/2)*(q<1)+((-(0.5/(q+0.001*(q==0)))*(2-q)**2))*((q>1)*(q<2)))
    return gradW
