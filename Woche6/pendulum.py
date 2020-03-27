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
    # TODO: Implement
    return 0.

def phiB(u, dt):
    # TODO: Implement
    return 0.

def energy(y):
    # TODO: Implement
    return 0.

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

# TODO: Implement

f = lambda tau,y: 0.

# RK45
t45, y45 = ode45(f, tspan, yo)

# Splitting Verfahren
methods = ['SS', 'PRKS6', 'Y61', 'KL8']
tSplit = []
ySplit = []

# TODO: Implement