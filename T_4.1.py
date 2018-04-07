#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Samuel Garcia

Solución del Problema 4.1 del Libro "An Introduction to Computational Fluid Dynamics"
de H K Versteeg and W Malalasekera
----------------------------------
El Planteamiento del Problema es el siguiente:
    
Considere el problema de la conducción de calor sin fuente en una varilla aislada 
cuyos extremos se mantienen a temperaturas constantes de 100 ° C y 500 ° C, respectivamente. 
El problema unidimensional esbozado en la Figura 4.3 se rige por
$
\dfrac{d}{dx} \left( k \dfrac{d T}{dx} \right) = 0
$

___________________________________________________________

 |<-------------- 0.5 m ------------->|
 A                                    B
 |====================================|
 |                                    |
 |====================================|
T_A = 100                            T_B = 500
___________________________________________________________
Figura 4.3

Calcule la distribución de temperatura en estado estable en la varilla. 
La conductividad térmica k es igual a 1000 W / m.K, el área de sección 
transversal A es 10 × 10-3 m2.

Solución Analitica:
    
T = 800 * x + 100
"""

#Importar Script donde se encuentra el Metodo de Volumen Finito
import FiniteVolumeMethod as fvm
#Librerias de Python
import numpy as np
import matplotlib.pyplot as plt
import os

# Datos del Problema (Variables)
longitud = 0.5 # metros
TA = 100 # °C 
TB = 500 # °C 
k  = 1000 # W/m.K
N  = 6 # Número de nodos

#
# Creamos la malla y obtenemos datos importantes
#
malla = fvm.Mesh(nodes = N, length = longitud)
nx    = malla.nodes()     # Número de nodos
nvx   = malla.volumes()   # Número de volúmenes
delta = malla.delta()     # Tamaño de los volúmenes

#
# Imprimimos los datos del problema
#
fvm.printData(Longitud = longitud,
              Temperatura_A = TA,
              Temperatura_B = TB,
              Conductividad = k,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)

#
#  Creamos los coeficientes del Metodo de Volumen Finito
#
df1 = fvm.Diffusion1D(nvx, Gamma = k, dx = delta) #Invocacion del Metodo de Volumen Finito
df1.alloc(nvx) # Se aloja memoria para los coeficientes
df1.calcCoef() # Se calculan los coeficientes

#
# Se construye el arreglo donde se guardará la solución
#
T = np.zeros(nvx) # El arreglo contiene ceros
T[0]  = TA        # Condición de frontera izquierda
T[-1] = TB        # Condición de frontera derecha
df1.bcDirichlet('LEFT_WALL', T[0])   # Se actualizan los coeficientes
df1.bcDirichlet('RIGHT_WALL', T[-1]) # de acuerdo a las cond. de frontera

# Se construye el sistema lineal de ecuaciones a partir de los coeficientes
# del Metodo de Volumen Finito
#
Su = df1.Su()  # Vector del lado derecho
A = fvm.Matrix(malla.volumes())  # Matriz del sistema
A.build(df1) # Construcción de la matriz en la memoria

#
# Se resuelve el sistema usando un algoritmo del módulo linalg
#
T[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
# Impresión del vector Solución de la Matriz construida
print('Solución = {}'.format(T))
print('.'+'-'*70+'.')
#
# Se construye un vector de coordenadas del dominio
#
x = malla.createMesh()
#
# Calculamos la solución analítica
#
Ta = 800 * x + 100
#
#  Se grafica la solución
#
x *= 100 # Transformación a [cm]
plt.plot(x,Ta, '-', label = 'Sol. analítica') # Sol. analítica
plt.plot(x,T,'o', label = 'Sol. FVM')
plt.title('Solución de $k (\partial^2 T/\partial x^2) = 0$ con FVM')
plt.xlabel('$x$ [cm]')
plt.ylabel('Temperatura [$^o$C]')
plt.grid()
plt.legend()
plt.savefig('Tarea_4_1.png')
plt.show()
os.system('Tarea_4_1.png')
os.system("pause")