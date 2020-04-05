#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 17:56:15 2020

@author: v1per
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

""" Kernaufgabe Runge-Kutta Verfahren 1 Serie 6 Tim Launer """

""" Wir müssen die Verfahren aus Teilaufgabe a) nun implementieren """

""" Wir haben ein Anfangswert gegeben und müssen dieses verwenden um die Konvergenzordnung
    der Methoden empirisch zu finden """
    
""" Anfangswerte """

y_start = np.array([[1],[0]])
T = 4.0
N = 100
    
""" Zielfunktion und exakte Lösung Implementieren """

def f(t,y):
    return y - 2*np.sin(t)

exact = lambda t: np.cos(t) + np.sin(t)
    
""" Step Methoden Implizit und explizit aus dem Skript """

def RK_expl_step(B, t, y, h, rhs, dim):
    """
    B: Butcher Schema
    t: time
    y: y(t-h)
    h:time step
    rhs: ODE
    dim: dimension of y. (we could calculate it directly in the code) """
     
    
    # number of increments s
    s = B.shape[0] - 1
    
    # three parts of the butcher scheme
    b = B[s][1:]; c = B[:,0][:-1]; A = B[:-1,1:]
    
    # "container" for the increments
    K = np.zeros((dim,s))
    
    #now we calculate the increments.
    #since we are using an explicit method, we just compute them in the order they appear in
    for i in range(s):
        K[:,i] = rhs(t + c[i]*h, y + h*np.dot(A[i,:].T, K.T))
    
    tnew = t + h
    ynew = y + h*np.squeeze(np.dot(K,b))
    
    return tnew, ynew

def RK_impl_step(B, t, y, h, rhs, dim):
    """
    B: Butcher Schema
    t: time
    y: y(t-h)
    h:time step
    rhs: ODE
    dim: dimension of y. (we could calculate it directly in the code) """
     
    
    # number of increments s
    s = B.shape[0] - 1
    
    # three parts of the butcher scheme
    b = B[s][1:]; c = B[:,0][:-1]; A = B[:-1,1:]
    
    # "container" for the increments
    K = np.zeros((dim,s))
    
    #now we calculate the increments.
    for i in range(s):
        
        K = K.ravel() #flatten the array for use in fsolve()
        def F(x):
            z = x.reshape(s); temp = np.zeros(s)
            for i in range(s):
                temp[i] = rhs(t + c[i]*h, y + h*np.dot(A[i,:].T,z.T))
            return temp.ravel()
            
        K = fsolve(lambda x: x - F(x), np.zeros(s).ravel())
        #finally we reshape K:
        K.shape = (s)
    
    tnew = t + h
    ynew = y + h*np.squeeze(np.dot(K,b))
    
    return tnew, ynew
    

""" Implemetiere die Methoden """

def expEuler(y0, T, N):
    """ y0 : Anfangswerte , T : Endzeit ,  N : Anzahl Schritte """
    
    t, h = np.linspace(0, T, N, retstep="True")
    y = np.zeros(N); y[0] = y0
    
    # Butcher Schema
    B  = np.array([[0.,0.],[0.,1.]])
    
    for k in range(N-1):
        y[k+1] = RK_expl_step(B, t[k], y[k], h, rhs = f, dim = 1)[1]
    
    return t,y

def impEuler(y0, T, N):
    """ y0 : Anfangswerte , T : Endzeit ,  N : Anzahl Schritte """
    
    t, h = np.linspace(0, T, N, retstep="True")
    y = np.zeros(N); y[0] = y0
    
    # Butcher Schema
    B  = np.array([[1.,1.],[1.,1.]])
    
    for k in range(N-1):
        y[k+1] = RK_impl_step(B, t[k], y[k], h, rhs = f, dim = 1)[1]
    
    return t,y


def iMP(y0, T, N):
    """ y0 : Anfangswerte , T : Endzeit ,  N : Anzahl Schritte """
    
    t, h = np.linspace(0, T, N, retstep="True")
    y = np.zeros(N); y[0] = y0
    
    # Butcher Schema
    B  = np.array([[0.5, 0.5],[1.,1.]])
    
    for k in range(N-1):
        y[k+1] = RK_impl_step(B, t[k], y[k], h, rhs = f, dim = 1)[1]
    
    return t,y



""" Funktionen ausführen """

if __name__ == "__main__":
    
    t, y = expEuler(y_start[0], T, N)
    t1, y1 = impEuler(y_start[0], T, N)
    t2, y2 = iMP(y_start[0], T, N)
    
    
    plt.figure()
    plt.plot(t, exact(t), 'r', label = "exact")
    plt.plot(t, y, 'b', label = "expEuler_RK")
    plt.plot(t, y1, 'g', label = "implEuler_RK")
    plt.plot(t, y2, 'y', label = "impl-MPR_RK")
    plt.legend()
    plt.grid()
    plt.ylabel("Lösungswerte")
    plt.xlabel("Zeit")
    plt.show()
    
    res = 2**np.arange(5, 12)
    errEE = np.empty(res.shape)
    errIMP = np.empty(res.shape)
    errIE = np.empty(res.shape)
    yEx = exact(T)
    
    for k, N_res in enumerate(res):
        tEE, yEE = expEuler(y_start[0], T, N_res)
        tIMP, yIMP = iMP(y_start[0], T, N_res)
        tIE, yIE = impEuler(y_start[0], T, N_res)
        
        errEE[k] = abs(yEE[-1] - yEx)
        errIMP[k] = abs(yIMP[-1] - yEx)
        errIE[k] = abs(yIE[-1] - yEx)
    
    plt.figure()
    plt.loglog(res, errEE, label = "EE")
    plt.loglog(res, errIMP, label = "iMP")
    plt.loglog(res, errIE, label = "IE")
    
    plt.loglog(res, (1.*res)**(-1), ':k', label = "$n^{-1}$")
    plt.loglog(res, (1.*res)**(-2), '-k', label = "$n^{-2}$")
    plt.xlabel("Anzahl Intervalle")
    plt.ylabel("$y(t) - y_N$")
    plt.grid()
    plt.legend()
    plt.show()
    

