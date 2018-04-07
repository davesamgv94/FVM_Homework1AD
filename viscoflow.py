"""
Created on Fri Mar 30 12:51:33 2018

VISCOFLOW : Visualization of Complex Flows

@author: Luis Miguel de la Cruz Salas

"""
#
# @Author : Luis M. de la Cruz Salas, 2018
#
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#plt.style.use('ggplot')

def grafica(x, phi, title = None, label = None, kind = None):
    plt.title(title)
    plt.xlabel('$x$ [m]')
    plt.ylabel('$\phi$ [...]')
    if kind:
        plt.plot(x,phi,kind,label=label,lw=2)
    else:
        plt.plot(x,phi,'--', label = label)

def show(filename = None):
    plt.legend()
    plt.grid()
    if filename:
        plt.savefig(filename)
    plt.show()

if __name__ == '__main__':

    from scipy import special
    x = np.linspace(-3,3)
    grafica(x,special.erfc(x),title='Función ERFC', label='erfc')
#	show('hola.pdf')
    show()
    
    fig, ax = plt.subplots()
    x = np.arange(0,2*np.pi,0.01)
    line, = ax.plot(x,np.sin(x))
    
    def animate(i):
        line.set_ydata(np.sin(x+i/10.0))
        return line
    ani = FuncAnimation(fig, animate, np.arange(1,10), interval=100)
    plt.show()
    