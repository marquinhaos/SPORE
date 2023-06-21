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
     (o)  |   PROP-LAMDA   | (\'/)
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
angulos_3=[15,30,45,60]
gasto_3=[4.285,5.187,6.291,7.746]
orden_aprox=3

inter = interp1d(angulos_3, gasto_3, kind='cubic')

x_graf=np.linspace(angulos_3[0],angulos_3[-1],1000) #mapeamos el esprectro de la funcion para verla entera
y_graf=inter(x_graf)

coeffs=np.polyfit(angulos_3,gasto_3,orden_aprox) #expresion analitica de la funcion
aprox=np.poly1d(coeffs)

# Gráfica de la función y los puntos
plt.figure()
plt.title("Gasto de propulsante frente al ángulo T=3", fontsize = 14, color = 'black')
plt.plot(angulos_3, gasto_3, 'ro', label='Datos')
plt.plot(x_graf, y_graf, 'b-', label='Interpolación')
plt.plot(x_graf, aprox(x_graf), 'g-', label=aprox)
plt.legend(loc='lower right')
plt.xlabel('angulos [º]')
plt.ylabel('gasto [kg]')


angulos_9=[30,45,60,75]
gasto_9=[7.236,7.627,8.207,9.141]
orden_aprox=3

inter = interp1d(angulos_9, gasto_9, kind='cubic')

x_graf=np.linspace(angulos_9[0],angulos_9[-1],1000) #mapeamos el esprectro de la funcion para verla entera
y_graf=inter(x_graf)

coeffs=np.polyfit(angulos_9,gasto_9,orden_aprox) #expresion analitica de la funcion
aprox=np.poly1d(coeffs)

# Gráfica de la función y los puntos
plt.figure()
plt.title("Gasto de propulsante frente al ángulo T=9", fontsize = 14, color = 'black')
plt.plot(angulos_9, gasto_9, 'ro', label='Datos')
plt.plot(x_graf, y_graf, 'b-', label='Interpolación')
plt.plot(x_graf, aprox(x_graf), 'g-', label=aprox)
plt.legend(loc='lower right')
plt.xlabel('angulos [º]')
plt.ylabel('gasto [kg]')

      
plt.show() 