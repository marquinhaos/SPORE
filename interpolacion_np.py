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
     (o)  | INTERPOLACION  | (\'/)
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

angulo_min=15 #no debe ser 0 ya que dará error
angulo_max=105 #dependiendo del paso puede no corresponder con los 90º
paso=1
orden_aprox=6 #orden del polinomio de la aproximacion


#-------------------------------------------------------------------------------------------------------------------------------------#
#SOLVER
#-------------------------------------------------------------------------------------------------------------------------------------#
t1=time.time()
angulos_grados=np.arange(angulo_min,angulo_max,paso) #no lo ponemos directamente en radianes porque los queremos en grados para la grafica
angulos=(np.pi/180)*angulos_grados #vector de angulos
tiempos=np.array([])
for i in range(np.size(angulos)):
    alpha=angulos[i]
    V=ci.u*np.array([np.cos(alpha)*np.sin(ci.beta),np.sin(alpha)*np.sin(ci.beta),np.cos(ci.beta)])
    Yin=np.concatenate((ci.R,V,[cuerpo.m]))

    sol = solve_ivp(Fnp, (ci.t0,ci.tf), Yin, t_eval=ci.t, method='DOP853', atol=ci.tola, rtol=ci.tolr)          
    Y = sol.y

    #calculo de la altitud
    h=np.sqrt(Y[0]*Y[0]+Y[1]*Y[1]+Y[2]*Y[2])

    for j in range(np.size(h)): 
        if h[j]>cuerpo.Rsoi and ci.t[j]>0:
            tiempos=np.append(tiempos,ci.t[j])
            break

inter = interp1d(angulos_grados, tiempos, kind='cubic')

x_graf=np.linspace(angulos_grados[0],angulos_grados[-1],1000) #mapeamos el esprectro de la funcion para verla entera
y_graf=inter(x_graf)

coeffs=np.polyfit(angulos_grados,tiempos,orden_aprox) #expresion analitica de la funcion
aprox=np.poly1d(coeffs)

# Gráfica de la función y los puntos
plt.title("Tiempo de salida de Rsoi frente al ángulo", fontsize = 14, color = 'black')
plt.plot(angulos_grados, tiempos, 'ro', label='Datos')
plt.plot(x_graf, y_graf, 'b-', label='Interpolación')
plt.plot(x_graf, aprox(x_graf), 'g-', label=aprox)
plt.legend(loc='best')
plt.xlabel('angulos [º]')
plt.ylabel('tiempo [s]')

t2 = time.time()        
print('tiempo de interpolacion [s]=',t2-t1)
print("Tiempo maximo de salido [s]: ", print(max(tiempos)))
print(tiempos)
      
plt.show() 
