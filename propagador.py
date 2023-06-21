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

for i in range(np.size(Y[6])): #revisa la solución para encontrar cuando m=cte -> se apaga el motor
    if Y[6,i]==Y[6,i+1] and ci.t[i]>ci.ton:
        for j in range(6):
            print("Valor de la variable ", labels[j], units[j], " al apagar propulsión:", Y[j,i])
        print("Tiempo al apagar propulsión [s]: ", ci.t[i])
        break


#calculo del gasto de propulsante, tecnicamente lo tenemos definido con la masa final, pero de esta manera es mas facil de leer
gp=1500-Y[6]
print("Gasto total de propulsante [kg]: ", gp[-1])

#calculo de la altitud
alt=np.sqrt(Y[0]*Y[0]+Y[1]*Y[1]+Y[2]*Y[2])-cuerpo.Ri
#-------------------------------------------------------------------------------------------------------------------------------------#
#SOLVER II
#-------------------------------------------------------------------------------------------------------------------------------------#       
t1 = time.time()
sol2 = solve_ivp(F, (ci.t0,ci.tf), ci.Yin, t_eval=ci.t, method='DOP853', atol=ci.tola*2, rtol=ci.tolr*2)          
Y2 = sol2.y

#calculo de la altitud
alt2=np.sqrt(Y2[0]*Y2[0]+Y2[1]*Y2[1]+Y2[2]*Y2[2])-cuerpo.Ri

Error=Y2-Y
Erroralt=alt2-alt

t2 = time.time()        
print('------- Tiempo de ejec. [s]  Solver II =',t2-t1, " ------")


#printeemos los errores maximos en cada variable
for i in range(7):
    print("Error maximo en la variable ", labels[i], units[i],": ", max(Error[i]))

#-------------------------------------------------------------------------------------------------------------------------------------#
#PLOTS 2D
#-------------------------------------------------------------------------------------------------------------------------------------#

#errores
for i in range(len(labels)):
    plt.figure()
    plt.plot(ci.t, Error[i,:], color = 'black', label=labels[i]) 
    plt.xlabel('time [s]')
    plt.ylabel(units[i])
    plt.title(r"Error en la variable "+labels[i], fontsize = 14, color = 'gray')
    plt.legend()

#variables
for i in range(len(labels)):
    plt.figure()
    plt.plot(ci.t, Y[i,:], color = 'black', label=labels[i]) 
    plt.xlabel('time [s]')
    plt.ylabel(units[i])
    plt.title(r"Solucion a la variable "+labels[i], fontsize = 14, color = 'gray')
    plt.legend()

#altitud
plt.figure()
plt.plot(ci.t, alt, color = 'black',label="altitud") 
plt.xlabel('time [s]')
plt.ylabel("km")
plt.title("Variacion de la altitud", fontsize = 14, color = 'gray')
plt.legend()

#error altitud
plt.figure()
plt.plot(ci.t, Erroralt, color = 'black',label="error actitud") 
plt.xlabel('time [s]')
plt.ylabel("km")
plt.title("error de la altitud", fontsize = 14, color = 'gray')
plt.legend()

#gasto de propulsante
plt.figure()
plt.plot(ci.t, gp, color = 'black',label="Gasto de propulsante") 
plt.xlabel('time [s]')
plt.ylabel("gasto [kg]")
plt.title("Gasto de propulsante", fontsize = 14, color = 'gray')
plt.legend()
#-------------------------------------------------------------------------------------------------------------------------------------#
#PLOTS 3D
#-------------------------------------------------------------------------------------------------------------------------------------#
#ORBITA
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(Y[0], Y[1], Y[2], color= "red" , cmap='Greens')
ax.set_xlim([-30000,30000]) # Reescalado de la animación
ax.set_ylim([-30000,30000])
ax.set_zlim([-30000,30000])

#creamos la Tierra
phi = np.linspace(0, np.pi, 100)
theta = np.linspace(0, 2 * np.pi, 100)
x = cuerpo.Ri * np.outer(np.sin(phi), np.cos(theta))
y = cuerpo.Ri * np.outer(np.sin(phi), np.sin(theta))
z = cuerpo.Ri * np.outer(np.cos(phi), np.ones_like(theta))
ax.plot_surface(x, y, z, rstride=4, cstride=4, color='blue')

#creamos la esfera de influencia
phi = np.linspace(0, np.pi, 100)
theta = np.linspace(0, 2 * np.pi, 100)
xsoi = cuerpo.Rsoi * np.outer(np.sin(phi), np.cos(theta))
ysoi = cuerpo.Rsoi * np.outer(np.sin(phi), np.sin(theta))
zsoi = cuerpo.Rsoi * np.outer(np.cos(phi), np.ones_like(theta))
ax.plot_surface(xsoi, ysoi, zsoi, rstride=4, cstride=4, color='red', alpha=0.06)

#ANIMACION
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim([-30000,30000]) # Reescalado de la animación
ax.set_ylim([-30000,30000])
ax.set_zlim([-30000,30000])
line, = ax.plot([], [], [], color='red') #color

#ploteamos la Tierra
ax.plot_surface(x, y, z, rstride=4, cstride=4, color='blue')
ax.plot_surface(xsoi, ysoi, zsoi, rstride=4, cstride=4, color='red', alpha=0.06)

def update(num, data, line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])
    return line

# Crea la animación
data = np.array([Y[0], Y[1], Y[2]])
line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1], 'r-')
ani = animation.FuncAnimation(fig, update, frames=ci.nstep, fargs=(data, line), blit=False, interval=1)

plt.show()