# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 21:08:02 2018

@author: Samuel Garcia
"""

import os 

line = '-' * 70
print('.'+ line + '.')
print('|{:^70}|'.format('Tarea 1'))
print('.'+ line + '.')
print('|{:^70}|'.format(' Autor: Luis M. De la Cruz, Adaptaciones: Samuel Garcia'))
print('.'+ line + '.')
print('|{:^70}|'.format('Menu de Programas'))
print('.'+ line + '.')
print('|{:^70}|'.format('1. Ejercicio 4.1'))
print('.'+ line + '.')
print('|{:^70}|'.format('2. Ejercicio 4.2'))
print('.'+ line + '.')
print('|{:^70}|'.format('3. Ejercicio 4.3'))
print('.'+ line + '.')
print('|{:^70}|'.format('4. Ejercicio 5.1'))
print('.'+ line + '.')
print('|{:^70}|'.format('5. Ejercicio 3 (Temporal)'))
print('.'+ line + '.')

programa = int(input("Escribe el numero de programa que deseas ejecutar: "))

if programa == 1:
    os.system('python T_4.1.py')
elif programa == 2:
    os.system('python T_4.2.py')
elif programa == 3:    
    os.system('python T_4.3.py')
elif programa == 4:
    os.system('python T_5.1.py')
elif programa == 5:
    os.system('python T_T_3.py')
else:
    os.system("exit")
