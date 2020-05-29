# -*- encoding: utf-8 -*-

import numpy as np
import numpy.linalg
from numpy.linalg import norm
from numpy import dot

import scipy
import scipy.linalg

# Dies hilft mit der Formatierung von Matrizen. Die Anzahl
# Nachkommastellen ist 3 und kann hier ----------------------.
# geändert werden.                                           v
np.set_printoptions(linewidth=200, formatter={'float': '{: 0.3f}'.format})

ew = np.array([100.0, 10.0, 12.0, 0.04, 0.234, 3.92, 72.0, 42.0, 77.0, 32.0])
n = ew.size
Q, _ = np.linalg.qr(np.random.random((n, n)))
A = np.dot(Q.transpose(), np.dot(np.diag(ew), Q))

# TODO Finden sie den grössten Eigenwert vom A mittels Potenzmethode.

def power_method(A, x, maxit):
    
    #arrays for storing the approximation of the maximal eigenvector at every iteration    
    EW_k = []; EWrayleigh_k = []
    x_k = [] #container for x at every iteration
    k = 1; x = x / norm(x)
    
    while(k <= maxit):
        x = dot(A,x)
        x = x / norm(x)
        k = k + 1
        EWrayleigh_k += [dot(x, (dot(A,x))) / dot(x,x)]
        EW_k += [norm(dot(A,x)) / norm(x)]
        x_k += [x]
        
    EV = x #eigenvector
    EW = EW_k[-1] #eigenvalue
    return EV, EW, EW_k, EWrayleigh_k, x_k

maxiter = 50
x0 = np.random.rand(A.shape[0])
ev, max_ew, ew_k, eewray_k, x_k = power_method(A, x0, maxiter)

print('--' * 30)
print('Der maximale Eigenwert zu A nach %s Iterationen ist: ' % (maxiter))
print(max_ew)
print()
    

# TODO Finden sie den kleinsten Eigenwert vom A mittels inverser Potenzmethode.

def invit(A, tol, maxiter):
    
    LUP = scipy.linalg.lu_factor(A, overwrite_a=True) #pivoted LU decomposition of A
    n = A.shape[0] #dimension of the vector space
    x = np.random.rand(n) #random starting vector
    x /= np.linalg.norm(x) #normalize x
    scipy.linalg.lu_solve(LUP, x, overwrite_b=True) #solve using LU decomposition, update x via overwrite
    
    lold = 0
    lmin = 1. / np.linalg.norm(x)
    x *= lmin #normalize x
    #now we repeat the iteration
    
    if tol != None: 
        
        while(abs(lmin-lold) > tol*lmin):
            lold = lmin
            scipy.linalg.lu_solve(LUP, x, overwrite_b=True)
            lmin = 1./np.linalg.norm(x)
            x *= lmin
    
    else:
        for i in range(maxiter):
            lold = lmin
            scipy.linalg.lu_solve(LUP, x, overwrite_b=True)
            lmin = 1./np.linalg.norm(x)
            x *= lmin   
        
    return lmin

tol =  1e-10
maxiter = 30

min_ew = invit(A, tol, maxiter)

print('--' * 30)
print('Der minimale Eigenwert zu A mit Toleranz %s ist: ' % (tol))
print(min_ew)
print()
    

# TODO Finden sie den Eigenwert am nächsten bei 42.



