#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 15:32:38 2018

@author: Samuel Garcia

Problema 3 de la tarea

"""

#Importar Script donde se encuentra el Metodo de Volumen Finito
import FiniteVolumeMethod as fvm
#Librerias de Python
import numpy as np
import viscoflow as vfl
from scipy import special
import os

### Calculo Mediante la Solución Analitica
def analyticSol(x, u, t, Gamma):
 	divisor = 2 * np.sqrt(Gamma * t)
 	sol = 0.5 * (special.erfc((x - u * t)/ divisor) + 
 		np.exp(u * x) * np.exp(-Gamma) * special.erfc((x + u * t)/divisor))
 	return sol

title_graf = '$\partial \phi / \partial t + \partial(p u \phi)/\partial x= \partial (\Gamma \partial\phi/\partial x)/\partial x$'
# Datos del Problema (Variables)
L = 2.5 # m
rho = 1.0 # kg/m^3
u = 1.0 # m/s
Gamma = 0.001 # kg / m.s
phiA = 1 #
phiB = 0 #
N = 50 # Número de nodos
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

exac = analyticSol(x, u, delta_t * steps, Gamma)
vfl.grafica(x,exac,label='Sol. Exac',kind='b-')

for i in range(1,steps+1):
#
# Iniciamos con coeficientes cero
#
    print('Time step = {}'.format(i * delta_t), sep = '\t')
    coef.cleanCoefficients()
#
#  Calculamos los coeficientes de FVM de la Difusión
#
    DE,DW,DP = dif.calcCoef()

#
#  Calculamos los coeficientes de FVM de la Advección
#
    adv.setU(u)
    """Elegir entre los diferentes metodos"""
    adv.calcCoef('Centradas', phiA, phiB)
    #adv.calcCoef('Upwind', phiA, phiB)
    #adv.calcCoef('Quick', phiA, phiB)
    #adv.calcCoef('UpwindII', phiA, phiB)
#
#  Calculamos los coeficientes de FVM de la parte temporal
#
    tem.calcCoef(phi)
#
# Se aplican las condiciones de frontera
#
    """Elegir entre las diferentes condiciones"""
    coef.bcDirichlet('LEFT_WALL', phiA)   # Se actualizan los coeficientes
    coef.bcDirichlet('RIGHT_WALL', phiB) # de acuerdo a las cond. de frontera
    #coef.bcDirichlet_QuicK(phiA, phiB, delta, gamma = Gamma, DE = DE, DW = DW, rho=rho, u=u)
    #coef.bcDirichlet_LUD(phiA, phiB)
    ####################################################################################

#
# Se construye el sistema lineal de ecuaciones a partir de los coef. de FVM
#
    Su = coef.Su()  # Vector del lado derecho
    A = fvm.Matrix(malla.volumes())  # Matriz del sistema
    A.build(coef) # Construcción de la matriz en la memoria

# Se resuelve el sistema usando un algoritmo del módulo linalg
#
    phi[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])    
#
# Usamos Viscoflow para graficar
# VISCOFLOW : Visualization of Complex Flows
#
    if (i % 100 == 0):
        etiqueta = 'Step = {}'.format(i*delta_t)
        vfl.grafica(x,phi,title=title_graf, label=etiqueta)

vfl.show('Tarea_T_3.png')
os.system('Tarea_T_3.png')
os.system("pause")
