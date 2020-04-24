#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 15:17:20 2020

@author: v1per
"""

import numpy as np
import sympy as sp
from sympy import *  # maybe a bug? - sympy cannot execute certain derivatives if you do not import the module this way
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

""" Aufgabe 2 Nichtlineares System lösen """

# Räpresentiere Punkte als Tupel im R^2 - p = (x,y)

""" Gleichungen implementieren """

x, y = sp.symbols('x y')

f = sp.exp(x*y) + x**2 + y - 6/5
g = x**2 + y**2 + x - 11/20

dfdx = sp.diff(f, x, 1)
dfdy = sp.diff(f, y, 1)

dgdx = sp.diff(g, x, 1)
dgdy = sp.diff(g, y, 1)

dfdx = sp.lambdify((x,y), dfdx, 'numpy')
dfdy = sp.lambdify((x,y), dfdy, 'numpy')

dgdx = sp.lambdify((x,y), dgdx, 'numpy')
dgdy = sp.lambdify((x,y), dgdy, 'numpy')

f1 = sp.lambdify((x,y), f, 'numpy')
g1 = sp.lambdify((x,y), g, 'numpy')

""" User F an einem Punkt p, das wir auf Nullstellen untersuchen wollen """
F = lambda p : np.array([f1(p[0], p[1]), g1(p[0],p[1])])

""" Die Ableitung von F als Jacobi Matrix DF an einem Punkt p"""
DF = lambda p : np.array([[dfdx(p[0],p[1]), dfdy(p[0],p[1])], [dgdx(p[0],p[1]), dgdy(p[0],p[1])]])



""" Teilaufgabe a) - Implementiere Newton Verfahren """


def newton(p0, F, DF, tol, N=1000):
    
    p = p0
    error = np.linalg.norm(F(p))
    i = 0
    approx_list = [p]
    
    while error > tol and i < N:
        
        fprime = DF(p)
        f = F(p)
        
        p = p - np.linalg.solve(fprime, f)
        error = np.linalg.norm(F(p))
        approx_list.append(p)
        
        i+=1
        
    return p, approx_list, i


""" Teilaufgabe b) Konvergenzordnung untersuchen """

def convergence(p0, F, DF, tol):
    res, x_list, niter = newton(p0, F, DF, tol)
    
    error = np.zeros((len(x_list), ))
    
    for k, x in enumerate(x_list):
        error[k] = np.linalg.norm(res - x)
        
    return x_list[1:], error[:-1]
        
def evaluate(p0, tol):
    
    p, error = convergence(p0, F, DF, tol)
    niter = len(p)
    
    print()
    print('Startwert p0: x = %.2f, y = %.2f' % (p0[0],p0[1]))
    print("%5s %8s %15s %17s %15s" % ('Iter.', 'X', 'Y', 'Error', 'Ordnung'))
    
    for i in range(niter):
        if i >= 1 and i < niter - 2:
            order = ((np.log(error[i+1]) - np.log(error[i]))
                     / (np.log(error[i]) - np.log(error[i-1])))
        else:
            order = 0.0

        pattern = "{:5d} {:15.7e} {:15.7e} {:15.7e} {:8.2f} \n"
        print(pattern.format(i+1, p[i][0], p[i][1], error[i], order))
    
    


if __name__ == '__main__':
    
    
    # STARTWERTE
    p0 = np.array([0.6, 0.5])
    tol = 1e-14
    
    x, x_list, niter = newton(p0, F, DF, tol)
    
    print('Teilaufgabe a)')
    print('Die NST: (x: %.5f, y: %.5f)' % (x[0],x[1]))
    print('Abbruch nach %s Iterationen\n' % (niter))
    print('Toleranz: ' + str(tol) + ', Verfahren ergab schlussendlich für F bei NST x: ')
    print(F(x))
    print()
    
    
    """ Create NST Plot """
    n = 500
    x = np.linspace(-1.5, 1.5, n)
    X, Y = np.meshgrid(x, x)
    Zf = f1(X, Y)
    Zg = g1(X, Y)
    z = np.array([[x[0],x[1]],[0,0]])
    z = np.zeros((500,500))
    
    # plt.figure()
    # CS = plt.contour(X, Y, Zf, [0], colors='g', label="NST von $f_1$")
    # plt.gca().clabel(CS, inline=1, fontsize=10)
    # plt.contour(X, Y, Zg, [0], colors='b', linestyle='--', label='NST von $g_1$')
    # plt.plot(x[0] , x[1], 'ro', label='gefundene NST')
    # plt.legend()
    # plt.grid()
    # plt.show()
    # # plt.pcolormesh(X, Y, Zf, vmax=20, vmin=-20)
    # # plt.colorbar()
    # # plt.show()
    
    print('Teilaufgabe b)')
    
    p0 = np.array([0.6,0.5])
    evaluate(p0, tol)
    print()
    p0 = np.array([0.4,0.25])
    evaluate(p0, tol)
    print()
    p0 = np.array([-4.6,8.2])
    evaluate(p0, tol)
    print()
    
    




