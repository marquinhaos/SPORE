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
          |                \  (o)
     (o)  |     ON-OFF     | (\'/)
    (\'/)/                 \;:;
     :;  |                 /)
      ;: `-.._    /__..--'\.' ;:
          :;  `--' :;   :;
"""

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import sys


from funciones import Fnp
import condiciones_iniciales as ci
import cuerpos as cuerpo


#-------------------------------------------------------------------------------------------------------------------------------------#
#CARACTERISTICAS DE LA INTERPOLACION
#-------------------------------------------------------------------------------------------------------------------------------------# 

angulos=[15,30,45,60,75]
angulos_3=[15,30,45,60]
angulos_9=[30,45,60,75]
tiempos_off=[268500,294348,272688,207192,111288]
tiempos_3=[158748,192144,233016,286908]
tiempos_9=[89352,94164,101340,112860]

# Gráfica de la función y los puntos
plt.figure()
plt.title("Tiempos de captura frente tiempo de paso libre", fontsize = 14, color = 'black')
plt.plot(angulos_3, tiempos_3, 'r', label='Tiempos-3')
plt.plot(angulos_9, tiempos_9, 'b', label='Tiempos-9')
plt.plot(angulos, tiempos_off, 'g', label="Tiempos-off")
plt.legend(loc='upper right')
plt.xlabel('Angulos [º]')
plt.ylabel('Tiempo [s]')

      
plt.show() 