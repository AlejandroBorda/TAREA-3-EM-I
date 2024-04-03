########################################
########################################
#Tarea 3 EM
#Solucion general computacional
########################################
########################################

import numpy as np
import scipy 
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pylab as plt
import time
from math import sinh

#INPUTS/CONSTANTES
N=int(input("Numero de filas: ")) # filas
M=int(input("Numero de columnas: ")) # columnas
dx=float(input("Longitud paso dx=dy: "))
dy=dx
Vx0=float(input("Condicion de frontera x=0: "))
Vxa=float(input("Condicion de frontera x=a: "))
Vy0=float(input("Condicion de frontera y=0: "))
Vyb=float(input("Condicion de frontera y=b: "))
a=N*dx
b=M*dy
av=(Vy0 + Vyb + Vx0 + Vxa)/4

def c_sol(N,M,dx,dy,Vx0,Vxa,Vy0,Vyb,a,b,av):
    """ c_sol: funcion de solucion computacional por metodo de relajacion
    Variables:
    N: numero de filas
    M: numero de columnas
    dx: longitud paso en x
    dy: longitud paso en y
    Vx0, Vxa, Vy0, Vyb: condiciones de frontera para x=0, x=a, y=0, y=b respetctivamente
    a: maximo del dominio en x 
    b: maximo del dominio en b
    av: valor promedio de las condiciones de frontera

    A continuacion se explican las diferentes partes del programa a medida que se implementan
    """

    """Se define el dominio discretizado con el valor promedio av en todos los puntos """
    grid = np.full((N, M), av)
    
    """Se aplican las condiciones de frontera a sus valores correspondientes del dominio discretizado"""
    for i in range(0,M):
        grid[0][i]=Vy0
        grid[N-1][i]=Vyb
        
    for i in range(0,N):
        grid[i][0]=Vx0
        grid[i][M-1]=Vxa

    """ Se define el valor de convergencia para el cual, si una iteracion del metodo de relajacion tiene
    un cambio menor a este valor, se finaliza el metodo y se terminan las iteraciones como consecuencia de
    la convergencia, ademas se define un numero maximo de iteraciones para el caso de ausencia de convergencia """
    convergence_threshold = 1e-4
    max_iterations = 500
    iteration = 0

    """ Implementacion del metodo de relajacion, en cada iteracion de este ciclo para cada elemento interno del 
    dominio se reemplaza el valor de esta posicion por el valor del promedio de las 4 posiciones adyacentes 
    guardadas en las variables Fp, Fm, Cp, Cm, el valor nuevo se guarda en la variable nv
    ademas se incluye un condicional para aquellos valores que posiblemente no tengan 4 valores adyacentes.
    Ademas, las iteraciones van en el rango 1->N-2 y 1->M-1 para solo incluir los valores interiores de la red.
    """
    while iteration < max_iterations:
        max_change = 0

        for i in range(1, N-1):
            for j in range(1, M-1):
                div_count = 0
                nv = 0
                Fp = 0
                Fm = 0
                Cp = 0
                Cm = 0
                if (j < M - 1):
                    Fp = grid[i][j+1]
                    div_count += 1
                if (j > 0):
                    Fm = grid[i][j-1]
                    div_count += 1
                if (i < N - 1):
                    Cp = grid[i+1][j]
                    div_count += 1
                if (i > 0):
                    Cm = grid[i-1][j]
                    div_count += 1

                previous_value = grid[i][j]
                nv = (Fp + Fm + Cp + Cm) / div_count
                grid[i][j] = nv
                change = abs(nv - previous_value)
                max_change = max(max_change, change)

        iteration += 1
        if max_change < convergence_threshold:
            print("Converged after", iteration, "iterations.")
            break
    else:
        print("Maximum number of iterations reached without convergence.")

    """Se grafica la red obtenida tras la aplicacion del metodo de relajacion"""
    plt.imshow(grid, origin='lower', cmap='viridis', extent=[0, b, 0, a])
    plt.colorbar(label='Potencial')
    plt.title('Solucion numerica')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()

        
c_sol(N,M,dx,dy,Vx0,Vxa,Vy0,Vyb,a,b,av)





