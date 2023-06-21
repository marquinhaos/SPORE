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
          |  CONDICIONES  \  (o)
     (o)  |    INICIALES   | (\'/)
    (\'/)/                 \;:;
     :;  |                 /)
      ;: `-.._    /__..--'\.' ;:
          :;  `--' :;   :;
"""

import cuerpos as cuerpo
import numpy as np
#-------------------------------------------------------------------------------------------------------------------------------------#
#DECLARAMOS CONDICIONES INICIALES
#-------------------------------------------------------------------------------------------------------------------------------------#
tola=2.23e-12
tolr=2.23e-12

c=20.
t0m=1.2e-7
T0=cuerpo.m*t0m*(3)  #0.00017999999999999998
uinf=0.1 #velocidad del satélite en el infinito[km/s]
u=np.sqrt(uinf*uinf+(2.*cuerpo.mu)/cuerpo.Rsoi)
alpha=15*np.pi/180
beta=90*np.pi/180
ton=47500

#Traducimos las CI a coordenadas cartesianas y lo almacenamos junto la masa en un vector
R=np.array([-cuerpo.Rsoi,0,0])
V=u*np.array([np.cos(alpha)*np.sin(beta),np.sin(alpha)*np.sin(beta),np.cos(beta)])
Yin=np.concatenate((R,V,[cuerpo.m]))
#-------------------------------------------------------------------------------------------------------------------------------------#
#CREAMOS ARRAY TEMPORAL
#-------------------------------------------------------------------------------------------------------------------------------------#
nstep = 100000
t0 = 0.
tf = 1200000.   
dt = (tf-t0)/nstep
t = np.linspace(t0,tf,nstep+1, endpoint = True) 