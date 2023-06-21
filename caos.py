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
     (o)  |      CAOS      | (\'/)
    (\'/)/                  \;:;
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


from funciones import F
import condiciones_iniciales as ci
import cuerpos as cuerpo

labels = ['x', 'y', 'z', 'vx','vy','vz',"m"]
units = ['[km]', '[km]', '[km]', '[km/s]','[km/s]','[km/s]',"[kg]"] #usaremos estes vectores para graficas y printeo de info
#-------------------------------------------------------------------------------------------------------------------------------------#
#SOLVER
#-------------------------------------------------------------------------------------------------------------------------------------#
t1 = time.time()        
sol = solve_ivp(F, (ci.t0,ci.tf), ci.Yin, t_eval=ci.t, method='DOP853', atol=ci.tola, rtol=ci.tolr)          
Y = sol.y
t2 = time.time()        
print('------ Tiempo de ejec. [s]  Solver I =',t2-t1, " ------")

#-------------------------------------------------------------------------------------------------------------------------------------#
#CI II
#-------------------------------------------------------------------------------------------------------------------------------------#   
Yin=ci.Yin+(0,0,0,0,0,0,-100)
#-------------------------------------------------------------------------------------------------------------------------------------#
#SOLVER II
#-------------------------------------------------------------------------------------------------------------------------------------#       
t1 = time.time()
sol2 = solve_ivp(F, (ci.t0,ci.tf), Yin, t_eval=ci.t, method='DOP853', atol=ci.tola, rtol=ci.tolr)          
Y2 = sol2.y
t2 = time.time()        
print('------- Tiempo de ejec. [s]  Solver II =',t2-t1, " ------")

#-------------------------------------------------------------------------------------------------------------------------------------#
#CI III
#-------------------------------------------------------------------------------------------------------------------------------------#   
Yin=ci.Yin+(1000,1000,0,0,0,0,0)
#-------------------------------------------------------------------------------------------------------------------------------------#
#SOLVER III
#-------------------------------------------------------------------------------------------------------------------------------------#       
t1 = time.time()
sol3 = solve_ivp(F, (ci.t0,ci.tf), Yin, t_eval=ci.t, method='DOP853', atol=ci.tola, rtol=ci.tolr)          
Y3 = sol3.y
t2 = time.time()        
print('------- Tiempo de ejec. [s]  Solver III =',t2-t1, " ------")

#-------------------------------------------------------------------------------------------------------------------------------------#
#CI IV
#-------------------------------------------------------------------------------------------------------------------------------------#   
Yin=ci.Yin+(500,500,0,0,0,0,0)
#-------------------------------------------------------------------------------------------------------------------------------------#
#SOLVER IV
#-------------------------------------------------------------------------------------------------------------------------------------#       
t1 = time.time()
sol4 = solve_ivp(F, (ci.t0,ci.tf), Yin, t_eval=ci.t, method='DOP853', atol=ci.tola, rtol=ci.tolr)          
Y4 = sol4.y
t2 = time.time()        
print('------- Tiempo de ejec. [s]  Solver IV =',t2-t1, " ------")

#-------------------------------------------------------------------------------------------------------------------------------------#
#PLOTS 2D
#-------------------------------------------------------------------------------------------------------------------------------------#

#variables
for i in range(len(labels)):
    plt.figure()
    plt.plot(ci.t, Y[i,:], color = 'black', label="base")
    plt.plot(ci.t, Y2[i,:], color = 'red', label="-100 kg")
    plt.plot(ci.t, Y3[i,:], color = 'blue', label="+1000 km")
    plt.plot(ci.t, Y4[i,:], color = 'green', label="+500 km")
    plt.ylabel(units[i])
    plt.title(r"Comparación ante perturbaciones de la variable "+labels[i], fontsize = 14, color = 'gray')
    plt.legend()

plt.show()