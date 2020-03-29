# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt

from ode45 import ode45
from simple_splitting import integrate, strang_splitting
from splitting_parameters import splitting_parameters

def splitting_step(Phi_a, Phi_b, y0, dt, a, b):
    """Allgemeines Splittingverfahren, mit Parameter `a`, `b`.

    Input:
        Phi_a: Lösungsoperator A, Signatur `(y0, dt)`.
        Phi_a: Lösungsoperator B, Signatur `(y0, dt)`.

           y0: Aktuelle Approximation der Lösung.
           dt: wie üblich.

            a: Relative Schrittweite für `Phi_a`.
            b: Relative Schrittweite für `Phi_b`.

    Output:
        y1: Approximative Lösung der ODE zur Zeit `t+dt`.
    """

    # Tip: defensiv programmieren:
    if a.size != b.size:
        raise Exception("Dimensions mismatch![%s, %s]".format(str(a), str(b)))

    y = y0
    for s in range(a.size):
        y = Phi_a(y, a[s]*dt)
        y = Phi_b(y, b[s]*dt)

    return y

def splitting_method(Phi_a, Phi_b, y0, t_end, n_steps, key):
    a, b = splitting_parameters(key)
    method = lambda y, dt: splitting_step(Phi_a, Phi_b, y, dt, a, b)

    return integrate(method, y0, t_end, n_steps)

def PRKS6(Phi_a, Phi_b, y0, t_end, n_steps):
    return splitting_method(Phi_a, Phi_b, y0, t_end, n_steps, "PRKS6")

def Y61(Phi_a, Phi_b, y0, t_end, n_steps):
    return splitting_method(Phi_a, Phi_b, y0, t_end, n_steps, "Y61")

def KL8(Phi_a, Phi_b, y0, t_end, n_steps):
    return splitting_method(Phi_a, Phi_b, y0, t_end, n_steps, "KL8")

# --- Problem setup

def phiA(u, dt):
    # TODO: Implement done
    q = u[0]
    p = u[1]
    t = u[2] + dt
    return np.array([q + p*dt, p, t])

def phiB(u, dt):
    # TODO: Implement done
    q = u[0]
    if mu == 0:
       p = u[1] + A*np.sin(omega*u[2]) - dt*(np.sin(u[0])) 
    else:
       p = u[1] + (1/mu)*(np.sin(u[0]) - A*np.sin(omega*u[2]))*np.exp(-mu*dt) - (1/mu)*(np.sin(u[0]) - A*np.sin(omega*u[2])) 
    
    t = u[2] + dt
    
    return np.array([q, p, t])

def energy(y):
    # TODO: Implement
    Ekin = (0.5 * y[:,1]**2)
    Epot =  -np.cos(y[:,0])

    E = Ekin + Epot 
    return E

# --- Modelparameter ----
T = 200
g = 1.
l = 1.
A = 1.
omega = 1.3
mu = 0.1
qo = np.pi/3
po = 0
to = 0
yo = np.array([qo, po, to])
tspan = [0, T]

# TODO: Implement done

f = lambda tau,y: np.array([y[1], np.sin(y[0]) + A*np.cos(omega*y[2]) - mu*y[1], 1])
                            
# RK45
t45, y45 = ode45(f, tspan, yo)

# Splitting Verfahren
methods = ['SS', 'PRKS6', 'Y61', 'KL8']
tSplit = []
ySplit = []

# TODO: Implement done
# Approximiere Pendel mit verschiedenen Verfahren
for m in methods:
    print("Approximiere Pendel durch " + m + " mit Mu ungleich null")
    t, y = splitting_method(phiA, phiB, yo, T, 30000, m)
    tSplit.append(t)
    ySplit.append(y)


plt.figure()
plt.plot(t45, y45[:,0], label="RK45")
plt.title("Mu ungleich null")
for t, y, m in zip(tSplit, ySplit, methods):
    plt.plot(t, y[:,0], label = m)
plt.xlabel("Time")
plt.ylabel("Auslenkung")
plt.legend(loc="best")
plt.grid()

plt.figure()
plt.plot(t45, energy(y45), label = "RK45")
for t, y, m in zip(tSplit, ySplit, methods):
    plt.plot(t, energy(y), label = m)
plt.xlabel("Time")
plt.ylabel("Energy, mu nicht null")
plt.legend(loc="best")
plt.grid()


plt.figure()
plt.plot(y45[:,0], y45[:,1], label = "RK45")
for t, y, m in zip(tSplit, ySplit, methods):
    plt.plot(y[:,0], y[:,1], label = m)
plt.xlabel("Time")
plt.ylabel("Trajektorie")
plt.title("Mu ist nicht null")
plt.legend(loc="best")
plt.grid()



# Setze Mu gleich null und plotte nochmals
mu = 0
tSplit2 = []
ySplit2 = []
t452, y452 = ode45(f, tspan, yo)
for m in methods:
    print("Approximiere Pendel durch " + m + " mit Mu gleich null")
    t, y = splitting_method(phiA, phiB, yo, T, 30000, m)
    tSplit2.append(t)
    ySplit2.append(y)
    

plt.figure()
plt.plot(t452, y452[:,0], label="RK45")
plt.title("Mu gleich null")
for t, y, m in zip(tSplit2, ySplit2, methods):
    plt.plot(t, y[:,0], label = m)
plt.xlabel("Time")
plt.ylabel("Auslenkung")
plt.legend(loc="best")
plt.grid()


plt.figure()
plt.plot(t452, energy(y452), label = "RK45")
for t, y, m in zip(tSplit2, ySplit2, methods):
    plt.plot(t, energy(y), label = m)
plt.xlabel("Time")
plt.ylabel("Energy, mu null")
plt.legend(loc="best")
plt.grid() 

plt.figure()
plt.plot(y452[:,0], y452[:,1], label = "RK45")
for t, y, m in zip(tSplit2, ySplit2, methods):
    plt.plot(y[:,0], y[:,1], label = m)
plt.xlabel("Time")
plt.ylabel("Trajektorie")
plt.title("Mu ist null")
plt.legend(loc="best")
plt.grid()


    
    
    
    

    
    
    
    
    
    