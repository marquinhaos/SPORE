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
          | SatPOrtoREina  \  (o)
     (o)  |    project     | (\'/)
    (\'/)/   >>>SPORE<<<   \;:;
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


from funciones import Fnp
import condiciones_iniciales as ci
import cuerpos as cuerpo
from orbital_elements import cartesian_to_orbital_elements

labels = ['e', 'omega', 'p', 'nu', 'rp']
units = ['[º]', '[º]', '[º]', '[º]', '[º]'] #usaremos estes vectores para graficas y printeo de info
#-------------------------------------------------------------------------------------------------------------------------------------#
#SOLVER
#-------------------------------------------------------------------------------------------------------------------------------------#
t1 = time.time()        
sol = solve_ivp(Fnp, (ci.t0,ci.tf), ci.Yin, t_eval=ci.t, method='DOP853', atol=ci.tola, rtol=ci.tolr)          
Y = sol.y
t2 = time.time()        
print('------ Tiempo de ejec. [s]  Solver I =',t2-t1, " ------")

kep=cartesian_to_orbital_elements(cuerpo.mu, Y)
Y=np.array([[kep.e],[kep.omega],[kep.p],[kep.nu],[kep.rp]])
#-------------------------------------------------------------------------------------------------------------------------------------#
#SOLVER II
#-------------------------------------------------------------------------------------------------------------------------------------#       
t1 = time.time()
sol2 = solve_ivp(Fnp, (ci.t0,ci.tf), ci.Yin, t_eval=ci.t, method='DOP853', atol=ci.tola*2, rtol=ci.tolr*2)          
Y2 = sol2.y
t2 = time.time()   

kep2=cartesian_to_orbital_elements(cuerpo.mu, Y2)
Y2=np.array([[kep2.e],[kep2.omega],[kep2.p],[kep2.nu],[kep2.rp]])
print('------- Tiempo de ejec. [s]  Solver II =',t2-t1, " ------")

Error=Y2-Y
#printeemos los errores maximos en cada variable
for i in range(len(labels)):
    print("Error maximo en la variable ", labels[i], units[i],": ", np.max(Error[i,:]))
#-------------------------------------------------------------------------------------------------------------------------------------#
#tiempo de paso
#-------------------------------------------------------------------------------------------------------------------------------------#
e=(Y[0])[-1][-1]; p=(Y[2])[-1][-1]; nu=(Y[3])[-1][-1]
h=np.sqrt(p*cuerpo.mu)
Hin=2.*(np.arctanh(((e-1)/(e+1))^0.5)*np.tan(nu/2))
Mhin=e*np.sinh(Hin)-Hin
tpaso=2*((p^2)/h)*(((e^2)-1)^(-3/2))*Mhin
print(tpaso)

#-------------------------------------------------------------------------------------------------------------------------------------#
#PLOTS 2D
#-------------------------------------------------------------------------------------------------------------------------------------#

#errores
for i in range(len(labels)):
    plt.figure()
    plt.plot(ci.t, max(Error[i]), color = 'black', label=labels[i]) 
    plt.xlabel('time [s]')
    plt.ylabel(units[i])
    plt.title(r"Error en la variable "+labels[i], fontsize = 14, color = 'gray')
    plt.legend()

#variables
for i in range(len(labels)):
    plt.figure()
    plt.plot(ci.t, max(Y[i]), color = 'black', label=labels[i]) 
    plt.xlabel('time [s]')
    plt.ylabel(units[i])
    plt.title(r"Solucion a la variable "+labels[i], fontsize = 14, color = 'gray')
    plt.legend()

plt.show()