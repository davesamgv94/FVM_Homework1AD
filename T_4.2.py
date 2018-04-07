#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Samuel Garcia

Solución del Problema 4.2 del Libro "An Introduction to Computational Fluid Dynamics"
de H K Versteeg and W Malalasekera

----------------------------------
Ahora hablamos de un problema que incluye fuentes distintas de las derivadas de las condiciones de contorno. 
La figura 4.6 muestra una gran placa de espesor L = 2 cm con conductividad térmica constante k = 0.5 W / m.K 
y generación de calor uniforme q = 1000 kW / m3. Las caras A y B están a temperaturas de 100 ° C y 200 ° C, respectivamente. 
Suponiendo que las dimensiones en las direcciones y y z son tan grandes que los gradientes de temperatura 
son significativos solo en la dirección x, calcule la distribución de temperatura en estado estable. 
Compare el resultado numérico con la solución analítica. La ecuación gobernante es:
$
\dfrac{d}{dx} \left( k \dfrac{d T}{dx} \right) + q = 0
$

___________________________________________________________

 |<-------------- 2.0 cm ------------>|
 A                                    B
 |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
 |                                    |
 |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

T_A = 100                            T_B = 200
___________________________________________________________
Figura 4.6

Solución Analitica:
    
$T = \left( \frac{T_B - T_A}{L} + \frac{q}{2k}(L - x) \right) x + T_A$

"""

#Importar Script donde se encuentra el Metodo de Volumen Finito
import FiniteVolumeMethod as fvm
#Librerias de Python
import numpy as np
import matplotlib.pyplot as plt
import os

### Calculo Mediante la Solución Analitica
def analyticSol(x):
    return ((TB - TA) / longitud + q * (longitud - x) / (2 * k) ) * x + TA

# Datos del Problema (Variables)
longitud = 0.02 # meters
TA = 100  # °C 
TB = 200  # °C 
k  = 0.5  # W/m.K
q  = 1e+06 # 1e+06 W/m^3 = 1000 kW/m^3 Fuente uniforme
N  = 6 # Número de nodos


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
              Temperatura_B = TB,
              Conductividad = k,
              Fuente = q,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)
#
#  Creamos los coeficientes de FVM
#
df1 = fvm.Diffusion1D(nvx, Gamma = k, dx = delta)
df1.alloc(nvx) # Se aloja memoria para los coeficientes
df1.calcCoef() # Se calculan los coeficientes
df1.setSu(q)  # Se agrega la fuente uniforme
#
# Se construye el arreglo donde se guardará la solución
#
T = np.zeros(nvx) # El arreglo contiene ceros
T[0]  = TA        # Condición de frontera izquierda
T[-1] = TB        # Condición de frontera derecha
df1.bcDirichlet('LEFT_WALL', T[0])   # Se actualizan los coeficientes
df1.bcDirichlet('RIGHT_WALL', T[-1]) # de acuerdo a las cond. de frontera
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
Ta = analyticSol(x1)
#
#  Se grafica la solución
#
x *= 100 # Transformación a [cm]
x1 *= 100
plt.plot(x1,Ta, '-', label = 'Sol. analítica') 
plt.plot(x,T,'o', label = 'Sol. FVM')
plt.title('Solución de $k (\partial^2 T/\partial x^2)+q = 0$ con FVM')
plt.ylim(75,300)
plt.xlabel('$x$ [cm]')
plt.ylabel('Temperatura [$^o$C]')
plt.grid()
plt.legend()
plt.savefig('Tarea_4_2.png')
plt.show()
os.system('Tarea_4_2.png')
os.system("pause")
