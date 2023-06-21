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
          |     Inter     \  (o)
     (o)  |     gasto     | (\'/)
    (\'/)/                 \;:;
     :;  |                 /)
      ;: `-.._    /__..--'\.' ;:
          :;  `--' :;   :;
"""
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy.integrate as integrate
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import sys

import condiciones_iniciales as ci
import cuerpos as cuerpo

#-------------------------------------------------------------------------------------------------------------------------------------#
#INTEGRADOR
#-------------------------------------------------------------------------------------------------------------------------------------#
#No usamos el de las funciones porque necesitamos cambiar ton en cada iteracion

def F(t, Y):
    FF=np.zeros_like(Y)

    r=np.sqrt(Y[0]*Y[0]+Y[1]*Y[1]+Y[2]*Y[2]) #Modulo del radio
    u=Y[3]*Y[3]+Y[4]*Y[4]+Y[5]*Y[5] #Modulo de la velocidad al cuadrado

    epsilon=(u)/2-cuerpo.mu/r

    if epsilon>=-cuerpo.mu/cuerpo.Rsoi and t>ton: #usamos una variable auxiliar prop para encender y apagar la propulsion
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

#-------------------------------------------------------------------------------------------------------------------------------------#
#CARACTERISTICAS DE LA INTERPOLACION
#-------------------------------------------------------------------------------------------------------------------------------------#
ton_min=0 #no debe ser 0 ya que dará error
ton_max=500000
paso=10000
#-------------------------------------------------------------------------------------------------------------------------------------#
#SOLVER
#-------------------------------------------------------------------------------------------------------------------------------------#
gasto=np.array([])
tiempos=np.array([])
for ton in range (ton_min,ton_max,paso):
  t1 = time.time()        
  sol = solve_ivp(F, (ci.t0,ci.tf), ci.Yin, t_eval=ci.t, method='BDF', atol=ci.tola, rtol=ci.tolr) #usamos BDF porque buscamos rapidez y no tanta precisión     
  Y = sol.y
  t2 = time.time()        
  print('------ Tiempo de ejec. [s]  Solver ton= [',ton,'] =',t2-t1, " ------")
  gp=1500-Y[6]
  print("Gasto total de propulsante [kg]: ", gp[-1])
  gasto=np.append(gasto,gp[-1])
  tiempos=np.append(tiempos,ton)
print(min(gasto), min(tiempos))
# Gráfica de la función y los puntos
plt.title("Gasto frente al tiempo de encendido", fontsize = 14, color = 'black')
plt.plot(tiempos, gasto, 'ro', label='Datos')
plt.xlabel('tiempo [s]')
plt.ylabel('gasto [kg]')

plt.show()
