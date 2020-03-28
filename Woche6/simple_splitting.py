# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from ode45 import ode45

# Wir wollen
# dy/dt = A(t) y + b*y
# mittels Splitting lösen. Dabei ist A eine Rotationsmatrix und b
# eine reelle Zahl.

# Die Lösung ist also eine Komposition von Rotation
# und Streckung von y.

b = -0.1
theta = 0.25*np.pi

def rhs(t, y):
    return np.dot(dRdt(t), np.dot(invR(t), y)) + b*y

def R(t):
    angle = theta*t
    A = np.array([[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle), np.cos(angle)]])
    return A

def invR(t):
    return R(-t)

def dRdt(t):
    angle = theta*t
    A = theta*np.array([[-np.sin(angle), -np.cos(angle)],
                        [np.cos(angle), -np.sin(angle)]])

    return A

def rotate(y, dt):
    A = R(dt)

    return np.dot(A, y)

def Phi_rot(y, dt):
    """Exakter Lösungsoperator der Rotationsterme in der ODE."""

    return rotate(y, dt)

def Phi_stretch(y, dt):
    """Exakter Lösungsoperator der Zerfallsterme in der ODE."""

    return np.array(y*np.exp(b*dt))

def integrate(method, y0, t_end, n_steps):
    """Führt wiedeholt Zeitschritte des Verfahrens `method` aus.

    Input:
        method : Einschrittverfahren mit Signatur `(y[i,:], t, dt)`.
    """

    t, dt = np.linspace(0, t_end, n_steps+1, retstep=True)
    y = np.empty((n_steps+1,) + y0.shape)
    y[0,...] = y0

    for i in range(n_steps):
        y[i+1,...] = method(y[i,...], dt)

    return t, y

def strang_splitting_step(Phi_a, Phi_b, y0, dt):
    # TODO implement
    # Beim Strang splitting müssen wir die Formel verwenden, um die Evolutionsoperatoren
    # Phi_a und Phi_b zu kombinieren: Phi = Phi_a(Phi_b(Phi_a(y0, dt/2), dt), dt/2)
    
    return Phi_a(Phi_b(Phi_a(y0, dt/2), dt), dt/2) # FIXME done 

def strang_splitting(Phi_a, Phi_b, y0, t_end, n_steps):
    method = lambda y, dt : strang_splitting_step(Phi_a, Phi_b, y, dt)

    return integrate(method, y0, t_end, n_steps)

if __name__ == "__main__":
    y0 = np.array([1.0, 0.0])
    t_end = 100.0
    n_steps = 1000

    t, y = strang_splitting(Phi_rot, Phi_stretch, y0, t_end, n_steps)

    t_ode45, y_ode45 = ode45(rhs, [0, t_end], y0)

    plt.plot(y[:,0], y[:,1],label = "Strang")
    plt.plot(y_ode45[:,0], y_ode45[:,1],label = "ode45")
    plt.legend(loc='best')
    # plt.savefig("spiral.pdf")
    plt.show()
