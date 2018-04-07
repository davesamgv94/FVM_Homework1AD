#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 15:32:38 2018

@author: Samuel Garcia

Solución del Problema 4.3 del Libro "An Introduction to Computational Fluid Dynamics"
de H K Versteeg and W Malalasekera

----------------------------------
En el ejemplo final trabajado de este capítulo discutimos el enfriamiento de una 
aleta circular por medio de transferencia de calor convectiva a lo largo de su longitud. 
La convección da lugar a un término de pérdida de calor dependiente de la temperatura 
o sumidero en la ecuación gobernante. En la Figura 4.9 se muestra una aleta cilíndrica 
con un área transversal uniforme A. La base está a una temperatura de 100 ° C (T_A) y 
el extremo está aislado. La aleta está expuesta a una temperatura ambiente de 20 ° C. 
La transferencia de calor dimensional en esta situación se rige por

$
\dfrac{d}{dx} \left( k A \dfrac{d T}{dx} \right) - hP(T - T∞) = 0
$

donde h es el coeficiente de transferencia de calor por convección, P el perímetro, 
k la conductividad térmica del material y $ T_ \ infty $ la temperatura ambiente. 
Calcule la distribución de temperatura a lo largo de la aleta y compare los 
resultados con la solución analítica dada por

$\frac{T - T∞}{T_A - T∞} = \frac{\cosh[n(L-x)]}{\cosh(nL)}$

Donde $n^2 = hP/(kA)$, $L$ es la longitud de la aleta y $x$ la distancia a lo largo de la aleta. 
Datos: $L = 1$ m, $hP/(kA) = 25/m^2$ (Note que $kA$ es constante).

La ecuación que rige en el ejemplo contiene un término sumidero, -hP (T - T∞), 
la pérdida de calor por convección, que es una función de la temperatura local T.

___________________________________________________________

 |<-------------- 1.0 m ------------->|
 A                                    
 |
 |
 |------------------------------------|
 |------------------------------------| Flujo igual a zero
 |        Temp. Ambiente (T∞)                           
 |                                    |
                                      |
T_A = 100                  \frac{\partial T}{\partial dx} = q = 0
___________________________________________________________
Figure 4.9

Cuando $kA$ = es constante, la ecuación gobernante se puede escribir como

$
\dfrac{d}{dx} \left( \dfrac{d T}{dx} \right) - n^2(T - T∞) = 0
$

"""

#Importar Script donde se encuentra el Metodo de Volumen Finito
import FiniteVolumeMethod as fvm
#Librerias de Python
import numpy as np
import matplotlib.pyplot as plt
import os

### Calculo Mediante la Solución Analitica
def analyticSol(x):
    return (TA - Tambiente) * np.cosh(n * (longitud - x)) / np.cosh(n * longitud) + Tambiente

# Datos del Problema (Variables)
longitud = 1.0 # metros
Tambiente = 20  # °C 
TA = 100  # °C 
n2 = 25 # /m^2 = hP/(kA)
n = 5
fluxB = 0 # Flujo igual a cero
N = 6 # Número de nodos

#
# Creamos la malla y obtenemos datos importantes
#
malla = fvm.Mesh(nodes = N, length = longitud)
nx    = malla.nodes()     # Número de nodos
nvx   = malla.volumes()   # Número de volúmenes
delta = malla.delta()     # Tamaño de los volúmenes
#
# Imprimimos los datos del problema (nicely)
#
fvm.printData(Longitud = longitud,
              Temperatura_A = TA,
              Flujo_B = fluxB,
              n2 = n2,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)
#
#  Creamos los coeficientes de FVM
#
df1 = fvm.Diffusion1D(nvx, Gamma = 1, dx = delta)
df1.alloc(nvx) # Se aloja memoria para los coeficientes
df1.calcCoef() # Se calculan los coeficientes
df1.setSu(n2 * Tambiente)  # Se agrega la fuente uniforme
df1.setSp(-n2)

#
# Se construye el arreglo donde se guardará la solución
#
T = np.zeros(nvx) # El arreglo contiene ceros
T[0]  = TA        # Condición de frontera izquierda
df1.bcDirichlet('LEFT_WALL', T[0])   # Se actualizan los coeficientes
df1.bcNeumman('RIGHT_WALL', fluxB) # de acuerdo a las cond. de frontera
#
# Se construye el sistema lineal de ecuaciones a partir de los coef. de FVM
#
Su = df1.Su()  # Vector del lado derecho
A = fvm.Matrix(malla.volumes())  # Matriz del sistema
A.build(df1) # Construcción de la matriz en la memoria

#
# Se resuelve el sistema usando un algoritmo del módulo linalg
#
T[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
T[-1] = T[-2] # Condición de frontera tipo Neumman
print('Solución = {}'.format(T))
print('.'+'-'*70+'.')
#
# Se construye un vector de coordenadas del dominio
#
x = malla.createMesh()
#
# Calculamos la solución exacta y el error
#
Ta = analyticSol(x)
error = fvm.calcError(T, Ta)
datos = {'x(m)': x,
         'T(x)': T,
         'Analytic': Ta,
         'Error': error}
fvm.printFrame(datos)
print('||Error|| = ', np.linalg.norm(error))
print('.'+ '-'*70 + '.')
#
# Calculamos la solución exacta en una malla más fina para graficar
#
x1 = np.linspace(0,longitud,100)
n = np.sqrt(n2)
Ta = analyticSol(x1)
#
#  Se grafica la solución
#
plt.plot(x1,Ta, '-', label = 'Sol. analítica') 
plt.plot(x,T,'o', label = 'Sol. FVM')
plt.title('Solución de $\partial^2 T/\partial x^2 - hP(T-T_\infty) = 0$ con FVM')
plt.ylim(10,110)
plt.xlabel('$x$ [m]')
plt.ylabel('Temperatura [$^o$C]')
plt.grid()
plt.legend()
plt.savefig('Tarea_4_3.png')
plt.show()
os.system('Tarea_4_3.png')
os.system("pause")