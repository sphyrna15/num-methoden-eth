#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 18:32:43 2020

@author: v1per
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from sympy import *


""" Aufgabe 2 """

""" reproduziere Plots aus Skirpt """

""" Broyden-Quasi-Newton Verfahren """

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
        
    return p, i  # optionally return approx_list



""" Broyden Verfahren aus Skript"""

from scipy.linalg import lu_solve, lu_factor
from numpy import dot, zeros

def fastbroyd(x0, F, J, tol=1e-12, maxit=20):
    
    x = x0.copy()
    DF = J(x)
    lup = lu_factor(DF) #LU decomposition of J
    k = 0; s = lu_solve(lup,F(x)) #linear system can be solved faster with LU deco
    x -= s; f = F(x); sn = dot(s,s) #step, sn is the squared norm of the correction
    
    #containers for storing s and sn:
    dx = zeros((maxit,len(x)))
    dxn = zeros(maxit)
    dx[k] = s; dxn[k] = sn
    k += 1; #k is the number of the iteration
    
    w = lu_solve(lup,f) #f = F(starting value) (see above)
    
    #now we perform the iteration
    
    f = F(x)
    
    while sn > tol and k < maxit:
        
        w = lu_solve(lup,f) #f = F(starting value) (see above)
        #now we update the correction for the k-th step
        #using the sherman morrison woodbury formel
        
        for r in range(1,k):
            w += dx[r]*( dot(dx[r-1],w) )/dxn[r-1]
            
        z = dot(s,w)
        s = (1+z/(sn-z))*w
        sn = dot(s,s)
        dx[k] = s
        dxn[k] = sn
        x -= s
        f = F(x)
        k+=1 #update x and iteration number k
    
    return x, k
    #return the final value and the numbers of iterations needed
    
    
def test(F, J, x0, exact, tol):
    
    methods = [newton, fastbroyd]
    labels = ["Newton", "Broyden"]
    iterations = [6, 11]
    
    plt.figure()
    
    for method, label, it in zip(methods, labels, iterations):
        
        k = np.zeros(it)
        x = np.zeros_like(x0)
        error = np.zeros_like(k)
        F_norm = np.zeros_like(k)
        
        for i in range(it):
            x, j = method(x0, F, J, tol, i+1)
            k[i] = j
            res = x - exact
            error[i] = np.linalg.norm(res)
            F_norm[i] = np.linalg.norm(F(x))
        
        plt.semilogy(k, error, '-*', label=label)
        plt.semilogy(k, F_norm, '.', label=label + " $\Vert F(x^k) \Vert$" ) 
    
    plt.grid()    
    plt.legend(loc='best')
              
    
    
    
    
    
def bsp_392():
    
    x, y = symbols('x y')

    f = x**2 - y**4
    g = x - y**3
    
    dfdx = diff(f, x, 1)
    dfdy = diff(f, y, 1)
    
    dgdx = diff(g, x, 1)
    dgdy = diff(g, y, 1)
    
    dfdx = lambdify((x,y), dfdx, 'numpy')
    dfdy = lambdify((x,y), dfdy, 'numpy')
    
    dgdx = lambdify((x,y), dgdx, 'numpy')
    dgdy = lambdify((x,y), dgdy, 'numpy')
    
    f1 = lambdify((x,y), f, 'numpy')
    g1 = lambdify((x,y), g, 'numpy')
    
    """ User F aus Bsp 3.9.2 an einem Punkt pen """
    F = lambda p : np.array([f1(p[0], p[1]), g1(p[0],p[1])])
    
    """ Jacobi Matrix für F """
    DF = lambda p : np.array([[dfdx(p[0],p[1]), dfdy(p[0],p[1])], [dgdx(p[0],p[1]), dgdy(p[0],p[1])]])
    
    x0 = np.array([0.7, 0.7])
    exact = np.array([1., 1.])
    tol = 1e-10
        
    test(F, DF, x0, exact, tol)
    
    
def bsp_395():
    
    n = 1000
    h = 2 / n
    tol = n*2.e-5
    
    one = np.ones(n)
    b = np.arange(1,n+1)
    a = 1 / (np.sqrt(dot(one,b) - 1)) * (b - one)
    a.reshape(n,1)
    A = np.identity(n)
    
    for i in range(n):
        for j in range(n):
            A[i][j] += a[i]*a[j]
    
    x0 = 2*one + np.arange(n)*h
    
    
    F = lambda x : np.dot(np.diag(x), np.dot(A, x)) - b
    DF = lambda x : np.dot(np.diag(x), A) + np.diag(np.dot(A, x))
    
    # F_395 = lambda x : diag(x).dot(A.dot(x)) - b
    # J_395 = lambda x : diag(x).dot(A) + diag(A.dot(x))
        
    exact = fsolve(F, x0)
    
    test(F, DF, x0, exact, tol)
    
    
    
if __name__ == '__main__':
    
    bsp_392()
    bsp_395()
    
    

 