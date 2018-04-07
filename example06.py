#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 15:32:38 2018

@author: luiggi

Problema 3 de la Tarea

"""

import FiniteVolumeMethod as fvm
import numpy as np
from scipy import special

def analyticSol(x, u, t, Gamma):
 	divisor = 2 * np.sqrt(Gamma * t)
 	sol = 0.5 * (special.erfc((x - u * t)/ divisor) + 
 		np.exp(u * x) * np.exp(-Gamma) * special.erfc((x + u * t)/divisor))
 	return sol

title_graf = '$\partial \phi / \partial t + \partial(p u \phi)/\partial x= \partial (\Gamma \partial\phi/\partial x)/\partial x$ con FVM'
L = 2.5 # m
rho = 1.0 # kg/m^3
u = 1.0 # m/s
Gamma = 0.001 # kg / m.s
phiA = 1 #
phiB = 0 #
N = 100 # Número de nodos
delta_t = 0.002 # Paso de tiempo
steps = 500
#
# Creamos la malla y obtenemos datos importantes
#
malla = fvm.Mesh(nodes = N, length = L)
nx    = malla.nodes()     # Número de nodos
nvx   = malla.volumes()   # Número de volúmenes
delta = malla.delta()     # Tamaño de los volúmenes
x = malla.createMesh()    # Vector de coordenadas del dominio
#
# Se construye el arreglo donde se guardará la solución
#
phi = np.zeros(nvx) # El arreglo contiene ceros
phi[0]  = phiA       # Condición de frontera izquierda
phi[-1] = phiB       # Condición de frontera derecha
#
# Imprimimos los datos del problema (nicely)
#
fvm.printData(Longitud = L,
              Densidad = rho,
              Velocidad = u,
              Coef_Diff = Gamma,
              Prop_0 = phiA,
              Prop_L = phiB,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)
#
# Se aloja memoria para los coeficientes
#
coef = fvm.Coefficients()
coef.alloc(nvx)
#
#  Definimos el tipo de coeficientes de la ecuación
#
dif = fvm.Diffusion1D(nvx, Gamma = Gamma, dx = delta)
adv = fvm.Advection1D(nvx, rho = rho, dx = delta)
tem = fvm.Temporal1D(nvx, rho = rho, dx = delta, dt = delta_t)
#
# Inicializamos algunos módulos para la animación
#
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
fig = plt.figure(figsize=(8,4))           # Figuras
ax = plt.axes(xlim=(0, 3), ylim=(0, 1.2)) # Ejes
line, = ax.plot(x, phi, '--', label='FVM')
label = ax.text(2.6, 0.5, 'Time = {:>8.5f}'.format(0),
                ha='center', va='center',
                fontsize=12)
#
# Graficamos la solución exacta para t = 1.0 (500 pasos)
#
exac = analyticSol(x, u, delta_t * steps, Gamma)
ax.plot(x,exac,'b-',label='Sol. Exac',lw=2)
plt.xlabel('$x$ [m]')
plt.ylabel('$\phi$ [...]')
ax.legend()
#
# Función que calcula la solución con FVM como función del tiempo
#
def implicitSolver(i):
    time_step = i * delta_t
    print('time step = {}'.format(time_step))
#
# Iniciamos con coeficientes cero
#
    coef.cleanCoefficients()
#
#  Calculamos los coeficientes de FVM de la Difusión
#
    dif.calcCoef()
    
    dE,dW,dP = dif.calcCoef()
#
#  Calculamos los coeficientes de FVM de la Advección
#
    adv.setU(u)    # Definimos la velocidad
    adv.calcCoef('Quick', phiA, phiB)
#
#  Calculamos los coeficientes de FVM de la parte temporal
#
    tem.calcCoef(phi)
#
# Se aplican las condiciones de frontera
#
    #coef.bcDirichlet('LEFT_WALL', phiA)   # Se actualizan los coeficientes
    #coef.bcDirichlet('RIGHT_WALL', phiB) # de acuerdo a las cond. de frontera
    coef.bcDirichlet_QuicK(phiA, phiB, delta, Gamma, dE, dW, dP, rho, u)
#
# Se construye el sistema lineal de ecuaciones a partir de los coef. de FVM
#
    Su = coef.Su()  # Vector del lado derecho
    A = fvm.Matrix(malla.volumes())  # Matriz del sistema
    A.build(coef) # Construcción de la matriz en la memoria
#
# Se resuelve el sistema usando un algoritmo del módulo linalg
#
    phi[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
#
# Actualizamos la gráfica
#
    line.set_ydata(phi) # cambia los datos en la dirección y
    label.set_text('Step = {:>8d} \n Time = {:>8.5f}'.format(i, time_step))
    ax.set_title(title_graf)
#
# Función que controla todo el cálculo y la animación
#
anim = FuncAnimation(fig,            # La figura
                     implicitSolver, # la función que cambia los datos
                     interval=1,     # Intervalo entre cuadros en milisegundos
                     frames=steps+1,   # Cuadros
                     repeat=False)   # Permite poner la animación en un ciclo


