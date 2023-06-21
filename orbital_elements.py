#!/usr/bin/env python
##
##==============================================================
## CURSO DE MECANICA ANALITICA Y ORBITAL - Universidade de Vigo
## Prof. Daniele Tommasini
##==============================================================
## Content:
##
## class cartesian_to_orbital_elements
## for obtaining the orbital elements from position and velocity vectors
##
## function orbital_elements_to_cartesian(self, mu, e, h, Omegadeg, ideg, omegadeg, nudeg)
## return x, y, z, vx, xy, vz
##
## 
# # Class cartesian_to_orbital_elements lo he escrito en forma de clase de Python
# # y lo he pensado para que de algun resultado incluso cuando
# # los elementos orbitales clasicos no estan definidos
# # La rutina esta hecha de manera que funciona tambien si cada componente
# # de la posicion y de la velocidad es un array de valores (un valor por cada t).
# # Ademas, da resultados incluso para valores singulares (sin i=0):
# # en tal caso, la linea de los nodos y Omega no tienen significado
# # (aunque el codigo retorna el valor cero para Omega), y omega, que sigue dando
# # la posicion del periapsis, se puede calcular a partir de una
# # referencia cualquiera (en la practica la rutina lo calcula desde el eje x).
# # Notar que la rutina funciona para cualquier valor de e,
# # aunque para e=0 la eleccion de $\omega$ es arbitraria
# # (la rutina la pone igual a cero en ese caso).
# #
# # function orbital_elements_to_cartesian es aplicacion directa de las formulas
# #
# # WARNING! Los valores de la inclinacion cerca de 0º o 180º solo se distinguen de cero dentro de 1.e-6 grados
# # Es el precio a pagar por usar los elementos orbitales clasicos, que son singulares por esos valores
# # Eso es debido a que cos(i) para i proxima a 0 o pi difiere de 1 solo por i**2/2 (en radianes!), por tanto la minima i
# # que se puede apreciar como distinto de cero (o de pi) tiene que ser tal que i**2/2 = 2.2e-16
## (es decir, el machine epsilon). Por tanto, delta i > 2e-8 rad, o en grados delta i > 1.e-6 grados. 
# # Todos los valores de la inclinacion que difieren de 0º o de 180º solo por 1.e-6 grados se confunden con 0º o con 180º.
# # Eso ocurre inevitablemente en la rutina orbital_elements_to_cartesian

# # 
# # En otro fichero (ejemplo_uso_orbital_elements) hay ejemplos de uso.


import numpy as np
from numpy import sin, pi, cos, arctan, sqrt, tan, array, arccos, zeros_like, arctanh, sinh

class cartesian_to_orbital_elements_noarray:
    def __init__(self, mu, Y):
        radtodegree = 180./pi
        self.Y = Y
        self.rvec = array([self.Y[0],self.Y[1],self.Y[2]])
        self.vvec = array([self.Y[3],self.Y[4],self.Y[5]])
        rdotr = self.rvec[0]*self.rvec[0] + self.rvec[1]*self.rvec[1] + self.rvec[2]*self.rvec[2]
        self.r = sqrt(rdotr)
        vdotv = self.vvec[0]*self.vvec[0] + self.vvec[1]*self.vvec[1] + self.vvec[2]*self.vvec[2]
        self.v = sqrt(vdotv)
        rdotv = self.rvec[0]*self.vvec[0] + self.rvec[1]*self.vvec[1] + self.rvec[2]*self.vvec[2]
        self.vr = rdotv/self.r
        self.hvec = array([self.rvec[1]*self.vvec[2] - self.rvec[2]*self.vvec[1],self.rvec[2]*self.vvec[0] - self.rvec[0]*self.vvec[2],self.rvec[0]*self.vvec[1] - self.rvec[1]*self.vvec[0]])
        hdoth = self.hvec[0]*self.hvec[0] + self.hvec[1]*self.hvec[1] + self.hvec[2]*self.hvec[2]
        self.h = sqrt(hdoth)
        self.Wvec = self.hvec/self.h
        self.evec = ((vdotv-mu/self.r)*self.rvec - rdotv*self.vvec)/mu
        edote = self.evec[0]*self.evec[0] + self.evec[1]*self.evec[1] + self.evec[2]*self.evec[2]
        self.e = sqrt(edote)  
        self.i = arccos(self.Wvec[2])*radtodegree
        self.p = hdoth/mu
        self.rp = self.p/(1.+self.e) 
        self.nvec = array([-self.hvec[1],self.hvec[0],0.*self.hvec[2]])
        NdotN = self.nvec[0]*self.nvec[0] + self.nvec[1]*self.nvec[1] + self.nvec[2]*self.nvec[2]
        self.Pvec = zeros_like(self.evec)
        self.Qvec = zeros_like(self.evec)
        if (self.i < 2.23e-16 or self.i > 180.*(1. - 2.23e-16)):
            self.Omega = 0.
            if self.e < 2.23e-16:
                cosnu = self.rvec[0]/self.r
                if (1.-cosnu < 2.23e-16): cosnu = 1.
                if (1.+cosnu < 2.23e-16): cosnu = -1.
                nurad = arccos(cosnu)
                if self.rvec[1]*cos(self.i) < 0.:
                    nurad = 2.*pi - nurad
                self.nu = nurad*radtodegree
                self.omega = 0.
            else:
                self.Pvec = self.evec/self.e
                self.Qvec = array([self.Wvec[1]*self.Pvec[2] - self.Wvec[2]*self.Pvec[1],self.Wvec[2]*self.Pvec[0] - self.Wvec[0]*self.Pvec[2],self.Wvec[0]*self.Pvec[1] - self.Wvec[1]*self.Pvec[0]])
                cosomega = self.Pvec[0]
                self.omega = arccos(cosomega)*radtodegree
                if self.Pvec[1]*cos(self.i) < 0.:
                   self.omega = 360. - self.omega
                cosnu = (self.Pvec[0]*self.rvec[0] + self.Pvec[1]*self.rvec[1] + self.Pvec[2]*self.rvec[2])/self.r
                if (1.-cosnu < 2.23e-16): cosnu = 1.
                if (1.+cosnu < 2.23e-16): cosnu = -1.
                nurad = arccos(cosnu)
                if rdotv < 0.:
                    nurad = 2.*pi - nurad
                self.nu = nurad*radtodegree
        else:
            self.nvec = self.nvec/sqrt(NdotN) 
            if self.e < 2.23e-16:
                self.Pvec[:] = self.nvec[:]
            else:
                self.Pvec[:] = self.evec[:]/self.e
                cosomega = self.nvec[0]*self.Pvec[0] + self.nvec[1]*self.Pvec[1] + self.nvec[2]*self.Pvec[2]
                if (1.-cosomega < 2.23e-16): cosomega = 1.
                if (1.+cosomega < 2.23e-16): cosomega = -1.
                self.omega = arccos(cosomega)*radtodegree
                if self.Pvec[2] < 0.:
                   self.omega = 360. - self.omega
                
            self.Qvec[:] = array([self.Wvec[1]*self.Pvec[2] - self.Wvec[2]*self.Pvec[1],self.Wvec[2]*self.Pvec[0] - self.Wvec[0]*self.Pvec[2],self.Wvec[0]*self.Pvec[1] - self.Wvec[1]*self.Pvec[0]])
            self.Omega = arccos(self.nvec[0])*radtodegree
            if self.nvec[1] < 0.:
                self.Omega = 360. - self.Omega
            cosnu = (self.Pvec[0]*self.rvec[0] + self.Pvec[1]*self.rvec[1] + self.Pvec[2]*self.rvec[2])/self.r
            if (1.-cosnu < 2.23e-16): cosnu = 1.
            if (1.+cosnu < 2.23e-16): cosnu = -1.
            nurad = arccos(cosnu)
            if rdotv < 0.:
                nurad = 2.*pi - nurad
            self.nu = nurad*radtodegree
            
        if self.e < 1.:
            self.a = self.p/(1.-self.e**2) 
            self.b = self.p/sqrt(1.-self.e**2) 
            self.T = 2.*pi*sqrt(self.a*self.a*self.a/mu)   
            self.EORH = 2.*arctan(sqrt((1.-self.e)/(1.+self.e))*tan(nurad/2.))    
            if self.nu > 180.:
                self.EORH = self.EORH + 2.*pi
            Mrad = self.EORH - self.e*sin(self.EORH)    
            self.M = Mrad*radtodegree   
        elif self.e == 1.:
            AAA = tan(self.nu/2.)
            self.M = AAA/2. + AAA*AAA*AAA/6.
        else:
            self.EORH = 2.*arctanh(sqrt((self.e-1.)/(1.+self.e))*tan(nurad/2.))   
            self.M = self.e*sinh(self.EORH) - self.EORH
            self.a = self.p/(self.e**2-1.) 
            self.b = self.p/sqrt(self.e**2-1.) 




class cartesian_to_orbital_elements:
    def __init__(self, mu, Y):
        radtodegree = 180./pi
        self.Y = Y
        self.rvec = array([self.Y[0],self.Y[1],self.Y[2]])
        self.vvec = array([self.Y[3],self.Y[4],self.Y[5]])
        rdotr = self.rvec[0]*self.rvec[0] + self.rvec[1]*self.rvec[1] + self.rvec[2]*self.rvec[2]
        self.r = sqrt(rdotr)
        vdotv = self.vvec[0]*self.vvec[0] + self.vvec[1]*self.vvec[1] + self.vvec[2]*self.vvec[2]
        rdotv = self.rvec[0]*self.vvec[0] + self.rvec[1]*self.vvec[1] + self.rvec[2]*self.vvec[2]
        self.vr = rdotv/self.r
        self.hvec = array([self.rvec[1]*self.vvec[2] - self.rvec[2]*self.vvec[1],self.rvec[2]*self.vvec[0] - self.rvec[0]*self.vvec[2],self.rvec[0]*self.vvec[1] - self.rvec[1]*self.vvec[0]])
        hdoth = self.hvec[0]*self.hvec[0] + self.hvec[1]*self.hvec[1] + self.hvec[2]*self.hvec[2]
        self.h = sqrt(hdoth)
        self.Wvec = self.hvec/self.h
        self.evec = ((vdotv-mu/self.r)*self.rvec - rdotv*self.vvec)/mu
        edote = self.evec[0]*self.evec[0] + self.evec[1]*self.evec[1] + self.evec[2]*self.evec[2]
        self.e = sqrt(edote)  
        self.i = arccos(self.Wvec[2])*radtodegree
        self.p = hdoth/mu
        self.rp = self.p/(1.+self.e) 
        self.nvec = array([-self.hvec[1],self.hvec[0],0.*self.hvec[2]])
        NdotN = self.nvec[0]*self.nvec[0] + self.nvec[1]*self.nvec[1] + self.nvec[2]*self.nvec[2]
        self.Pvec = zeros_like(self.evec)
        self.Omega = zeros_like(self.i)
        self.omega = zeros_like(self.i)
        self.Qvec = zeros_like(self.evec)
        self.nu = zeros_like(self.i)
        self.M = zeros_like(self.i)   
        self.EORH = zeros_like(self.i)  
        self.a = zeros_like(self.i)   
        self.b = zeros_like(self.i)   
        self.T = zeros_like(self.i)   
        for k in range(self.e.shape[0]):
            if (self.i[k] < 2.23e-16 or self.i[k] > 180.*(1. - 2.23e-16)):
                if self.e[k] < 2.23e-16:
                    cosnu = self.rvec[0,k]/self.r[k]
                    if (1.-cosnu < 2.23e-16): cosnu = 1.
                    if (1.+cosnu < 2.23e-16): cosnu = -1.
                    nurad = arccos(cosnu)
                    if self.rvec[1,k]*cos(self.i[k]) < 0.:
                        nurad = 2.*pi - nurad
                    self.nu[k] = nurad*radtodegree
                else:
                    self.Pvec[:,k] = self.evec[:,k]/self.e[k]
                    self.Qvec[:,k] = array([self.Wvec[1,k]*self.Pvec[2,k] - self.Wvec[2,k]*self.Pvec[1,k],self.Wvec[2,k]*self.Pvec[0,k] - self.Wvec[0,k]*self.Pvec[2,k],self.Wvec[0,k]*self.Pvec[1,k] - self.Wvec[1,k]*self.Pvec[0,k]])
                    cosomega = self.Pvec[0,k]
                    self.omega[k] = arccos(cosomega)*radtodegree
                    if self.Pvec[1,k]*cos(self.i[k]) < 0.:
                       self.omega[k] = 360. - self.omega[k]
                    cosnu = (self.Pvec[0,k]*self.rvec[0,k] + self.Pvec[1,k]*self.rvec[1,k] + self.Pvec[2,k]*self.rvec[2,k])/self.r[k]
                    if (1.-cosnu < 2.23e-16): cosnu = 1.
                    if (1.+cosnu < 2.23e-16): cosnu = -1.
                    nurad = arccos(cosnu)
                    if rdotv[k] < 0.:
                        nurad = 2.*pi - nurad
                    self.nu[k] = nurad*radtodegree
            else:
                self.nvec[:,k] = self.nvec[:,k]/sqrt(NdotN[k]) 
                if self.e[k] < 2.23e-16:
                    self.Pvec[:,k] = self.nvec[:,k]
                else:
                    self.Pvec[:,k] = self.evec[:,k]/self.e[k]
                    cosomega = self.nvec[0,k]*self.Pvec[0,k] + self.nvec[1,k]*self.Pvec[1,k] + self.nvec[2,k]*self.Pvec[2,k]
                    self.omega[k] = arccos(cosomega)*radtodegree
                    if self.Pvec[2,k] < 0.:
                       self.omega[k] = 360. - self.omega[k]
                    
                self.Qvec[:,k] = array([self.Wvec[1,k]*self.Pvec[2,k] - self.Wvec[2,k]*self.Pvec[1,k],self.Wvec[2,k]*self.Pvec[0,k] - self.Wvec[0,k]*self.Pvec[2,k],self.Wvec[0,k]*self.Pvec[1,k] - self.Wvec[1,k]*self.Pvec[0,k]])
                self.Omega[k] = arccos(self.nvec[0,k])*radtodegree
                if self.nvec[1,k] < 0.:
                    self.Omega[k] = 360. - self.Omega[k]
                cosnu = (self.Pvec[0,k]*self.rvec[0,k] + self.Pvec[1,k]*self.rvec[1,k] + self.Pvec[2,k]*self.rvec[2,k])/self.r[k]
                if (1.-cosnu < 2.23e-16): cosnu = 1.
                if (1.+cosnu < 2.23e-16): cosnu = -1.
                nurad = arccos(cosnu)
                if rdotv[k] < 0.:
                    nurad = 2.*pi - nurad
                self.nu[k] = nurad*radtodegree
                
            if self.e[k] < 1.:
                self.a[k] = self.p[k]/(1.-self.e[k]**2) 
                self.b[k] = self.p[k]/sqrt(1.-self.e[k]**2) 
                self.T[k] = 2.*pi*sqrt(self.a[k]*self.a[k]*self.a[k]/mu) 
                self.EORH[k] = 2.*arctan(sqrt((1.-self.e[k])/(1.+self.e[k]))*tan(nurad/2.))  
                if self.nu[k] > 180.:
                    self.EORH[k] = self.EORH[k] + 2.*pi
                Mrad = self.EORH[k] - self.e[k]*sin(self.EORH[k])  
                self.M[k] = Mrad*radtodegree   
            elif self.e[k] == 1.:
                AAA = tan(self.nu[k]/2.)
                self.M[k] = AAA/2. + AAA*AAA*AAA/6.
            else:
                self.a[k] = self.p[k]/(self.e[k]**2-1.) 
                self.b[k] = self.p[k]/sqrt(self.e[k]**2-1.) 
                self.EORH[k] = 2.*arctanh(sqrt((self.e[k]-1.)/(1.+self.e[k]))*tan(nurad/2.))   
                self.M[k] = self.e[k]*sinh(self.EORH[k]) - self.EORH[k]
                

def orbital_elements_to_cartesian(mu, e, h, Omegadeg, ideg, omegadeg, nudeg):
        degreetorad = pi/180.
        i = ideg*degreetorad
        Omega = Omegadeg*degreetorad
        omega = omegadeg*degreetorad
        nu = nudeg*degreetorad
        p = h*h/mu
        u = omega + nu
        r = p/(1.+e*cos(nu))
        x = r*(cos(u)*cos(Omega) - sin(u)*cos(i)*sin(Omega))
        y = r*(cos(u)*sin(Omega) + sin(u)*cos(i)*cos(Omega))
        z = r*(sin(u)*sin(i))
        rdot = mu*e*sin(nu)/h
        xdot = rdot*x/r - (sin(u)*cos(Omega) + cos(u)*cos(i)*sin(Omega))*h/r
        ydot = rdot*y/r + (-sin(u)*sin(Omega) + cos(u)*cos(i)*cos(Omega))*h/r
        zdot = rdot*z/r + (cos(u)*sin(i))*h/r
        return x, y, z, xdot, ydot, zdot


############### ejemplo de lo que se puede hacer con una class:

"""
class test:
    def __init__(self, Y):
        self.Y = Y
        self.portres = 3.*Y
    def cuadrado(self):
        return (self.Y)**2
    def __call__(self, x):
        return self.cuadrado()*x

#### ejemplo de uso:
YY = array([2.,3.])
print(test(YY).portres)  # llamo algo calculado en el init para toda la clase (por tanto llamado con self.)
print(test(YY).cuadrado()) # llamo algo definido como funcion dentro de la clase
print(test(YY)(10.))  # llamo la funcion __call__ de la clase
### alternativamente, puedo llamar las tres operaciones anteriores definiendo primero tes = test(YY):
# tes = test(YY)       
# print(tes.portres)  # llamo algo calculado en el init para toda la clase (por tanto llamado con self.)
# print(tes.cuadrado()) # llamo algo definido como funcion dentro de la clase
# print(tes(10.))  # llamo la funcion __call__ de la clase
###############
"""