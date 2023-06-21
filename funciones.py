""""

              _......_
           .-'      _.-`.
         .'    .-. '-._.-`.
        / /\   |  \        `-..
       / / |   `-.'      .-.   `-.
      /  `-'            (   `.    `.
     |           /\      `-._/      \
     '    .'\   /  `.           _.-'|
    /    /  /   \_.-'        _.':;:/
  .'     \_/             _.-':;_.-'
 /   .-.             _.-' \;.-'
/   (   \       _..-'     |
\    `._/  _..-'          |
 `-.....-'/               |
          |   FUNCIONES   \  (o)
     (o)  |                | (\'/)
    (\'/)/                 \;:;
     :;  |                 /)
      ;: `-.._    /__..--'\.' ;:
          :;  `--' :;   :;
"""
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy.integrate as integrate
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import sys

import numpy as np
import cuerpos as cuerpo
import condiciones_iniciales as ci
import __main__ as m

#integrador
def F(t, Y):
    FF=np.zeros_like(Y)

    r=np.sqrt(Y[0]*Y[0]+Y[1]*Y[1]+Y[2]*Y[2]) #Modulo del radio
    u=Y[3]*Y[3]+Y[4]*Y[4]+Y[5]*Y[5] #Modulo de la velocidad al cuadrado

    epsilon=(u)/2-cuerpo.mu/r

    if epsilon>=-cuerpo.mu/cuerpo.Rsoi and t>ci.ton: #usamos una variable auxiliar prop para encender y apagar la propulsion
        #T=[0.,0.,0.] #descomentar para anular prop
        T=-ci.T0/np.sqrt(u)*np.array([Y[3],Y[4],Y[5]]) #vector de fuerza de empuje
        FF[6]=-ci.T0/ci.c
    else:
        T=[0.,0.,0.]
        FF[6]=FF[6]
        

    r3=(r*r*r)

    FF[0]=Y[3]
    FF[1]=Y[4]
    FF[2]=Y[5]

    FF[3]=(-cuerpo.mu*Y[0])/r3+T[0]/Y[6]
    FF[4]=(-cuerpo.mu*Y[1])/r3+T[1]/Y[6]
    FF[5]=(-cuerpo.mu*Y[2])/r3+T[2]/Y[6]

    if r>cuerpo.Rsoi:
        raise ValueError("Fallo en la captura en tiempo:",t)
    if r<cuerpo.Ri:
        raise ValueError("Colision con el cuerpo en tiempo:",t)
        
    

    return FF


#integrador sin propulsion
def Fnp(t, Y):

    FF=np.zeros_like(Y)
    r=np.sqrt(Y[0]*Y[0]+Y[1]*Y[1]+Y[2]*Y[2]) #Modulo del radio

    r3=(r*r*r)

    FF[0]=Y[3]
    FF[1]=Y[4]
    FF[2]=Y[5]
    FF[3]=(-cuerpo.mu*Y[0])/r3
    FF[4]=(-cuerpo.mu*Y[1])/r3
    FF[5]=(-cuerpo.mu*Y[2])/r3
    FF[6]=FF[6]

    return FF
    

''''
#integrador
def Fp(t, Y):
    FF=np.zeros_like(Y)

    r=np.sqrt(Y[0]*Y[0]+Y[1]*Y[1]+Y[2]*Y[2]) #Modulo del radio
    u=Y[3]*Y[3]+Y[4]*Y[4]+Y[5]*Y[5] #Modulo de la velocidad al cuadrado

    epsilon=(u)/2-cuerpo.mu/r

    if epsilon>=-cuerpo.mu/cuerpo.Rsoi: #usamos una variable auxiliar prop para encender y apagar la propulsion
        #T=[0.,0.,0.] #descomentar para anular prop
        T=-ci.T0/np.sqrt(u)*np.array([Y[3],Y[4],Y[5]]) #vector de fuerza de empuje
        FF[6]=-ci.T0/ci.c
    else:
        T=[0.,0.,0.]
        FF[6]=FF[6]
        

    r3=(r*r*r)

    FF[0]=Y[3]
    FF[1]=Y[4]
    FF[2]=Y[5]

    FF[3]=(-cuerpo.mu*Y[0])/r3+T[0]/Y[6]
    FF[4]=(-cuerpo.mu*Y[1])/r3+T[1]/Y[6]
    FF[5]=(-cuerpo.mu*Y[2])/r3+T[2]/Y[6]
    
    return FF
    '''