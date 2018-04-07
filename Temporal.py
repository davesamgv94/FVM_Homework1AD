#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 18:46:43 2018

@author: Samuel Garcia
"""
#Librerias de Python
import numpy as np
# Importar Script donde se encuentra la parte
# de los Coeficienetes
from Coefficients import Coefficients

class Temporal1D(Coefficients):
    
    """
    Esta clase realiza el cambio en los coefiecientes con el paso del tiempo ya que
    estos van cambiando, estos se tiene que ir actualizando uno tras otro, dependiendo
    del tiempo que se asigne.
    """  
    # Declaración de variables
    def __init__(self, nvx = None, rho = None, dx = None, dt = None):
        super().__init__(nvx)
        self.__nvx = nvx
        self.__rho = rho
        self.__dx = dx
        self.__dt = dt

    def __del__(self):
        del(self.__nvx)
        del(self.__rho)
        del(self.__dx)
        del(self.__dt)
    
    def deltaT(self):
        return self.__dt
    
    #Calculo de los Coeficientes
    def calcCoef(self, phi_old):
        aP = self.aP()
        Su = self.Su()
        rho = self.__rho
        dx_dt = self.__dx / self.__dt
        
        #Aqui ocurre la actualizacion o cambio de los coeficientes
        #dependiendo del tiempo que se asigne
        for i in range(1,self.__nvx-1):
            aP[i] += rho * dx_dt 
            Su[i] += phi_old[i] * dx_dt

# Parte Principal del Script donde se invica la clase e Impresión de Pruebas   
if __name__ == '__main__':
    
    nx = 6
    phi_old = np.sin(np.linspace(0,1,nx))
    print('-' * 20)  
    print(phi_old)
    print('-' * 20)  

    tf1 = Temporal1D(6, 1, 1, 1)
    tf1.alloc(6)
    tf1.calcCoef(phi_old)
    print(tf1.aP())
    print('-' * 20)  




