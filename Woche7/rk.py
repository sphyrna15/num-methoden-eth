# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import fsolve

class RungeKutta(object):
    """Basisklasse für alle Runge-Kutta Verfahren."""

    def __init__(self, A, b):
        """Runge-Kutta Löser mit dem Butcher Tableau (A, b).

        Wir implementieren nur konsistente Runge-Kutta Verfahren,
        d.h. solche mit

            c[i] = sum_j A[i,j]

        Input:
            A: Butcher Matrix A.
            b: Butcher Vektor b.
        """
        self.A = A.copy()
        self.b = b.copy()
        self.c = np.sum(A, axis=1)

        self.n_stages = A.shape[0]

    # Die method `__call__` erlaubt uns das object wie eine Funktion zu
    # benutzen:
    #   rk = RungeKutta(A, b)
    #   t, y = rk(f, y0, t_end, n_steps)
    def __call__(self, f, y0, t_end, n_steps):
        """ Berechnet die Approximation der Lösung.

        Input:
                  f : Rechte Seite der ODE.
                 y0 : Anfangsbedingung.
              t_end : Endzeit.
            n_steps : Anzahl Zeitschritte.
        """
        y = np.empty((n_steps+1, y0.size))
        t, dt = np.linspace(0, t_end, n_steps+1, retstep=True)

        y[0,:] = y0
        for i in range(n_steps):
            y[i+1,:] = self.step(f, y[i,:], t[i], dt)

        return t, y

    # NOTE Wir wollen drei Fälle unterscheiden: explizites RK, semi-implizites
    #   RK und voll-implizites RK. Die drei Fälle müssen unterschieden werden,
    #   da man bei explizitem RK nicht unnötig ein Gleichungssystem lösen will.
    #   Diese Funktionalität wird also in drei seperaten Unterklassen
    #   implementiert.
    # def step(self, t, y0, dt):
    #     # ...
    #     return y1


class ExplicitRungeKutta(RungeKutta):
    """Implementiert explizite Runge-Kutta Verfahren."""
    # NOTE Den Konstruktor brauchen wir nicht nochmals zu implementieren. Der
    #   default ist der Konstruktor von RungeKutta (also der Basisklasse).
    #   Dies gilt auch für alle anderen Funktionen.
    # def __init__(self, A, b, c):
    #     self.A = A.copy()
    #     # ...

    # NOTE auch __call__ gibt es schon, siehe oben.
    # def __call__(self, f, y0, t_end, n_steps):
    #     # ...
    #     return t, y

    def step(self, rhs, y0, t, dt):
        """ Führt einen expliziten RK Schritt aus.

        Input:
            rhs : Rechte Seite der ODE.
             y0 : Approximative Lösung zur Zeit `t`.
              t : Aktuelle Simulationszeit.
             dt : Zeitschritt, \Delta t.
        """
        # Diese Funktion berechnet einen Schritt y^{n} --> y^{n+1}.
        # Aber Namen wie (yn, ynp1) sind etwas kryptisch, deswegen (y0, y1), denn
        # eigentlich können wir das Problem als eine ODE mit Anfangsbedinung
        # `y0 = y0`, `t0 = t`, `t_end = t + dt` betrachten. Es besteht also keine
        # Verwirrungsgefahr :)
        #
        # Eventuell ist (y_current, y_next) eine Alternative, aber viel länger
        # auch verglichen mit `t`, `dt`.

        print("ExRK step t = {:.3e} ...".format(t), end="")

        A, b, c = self.A, self.b, self.c
        s = self.n_stages
        dydt = np.empty((y0.size, s))   # werden auf Papier `k` genannt.

        for i in range(s):
            y_star = y0 + dt*np.sum(A[i,:i]*dydt[:,:i], axis=-1)
            dydt[:,i] = rhs(t + dt*c[i], y_star)

        y1 = y0 + dt*np.sum(b*dydt, axis=-1)

        print(" done.")
        return y1


class SemiImplicitRungeKutta(RungeKutta):
    """Implementiert diagonal implizite Runge-Kutta Verfahren."""

    def step(self, rhs, y0, t, dt):
        """Berechnet einen Schritt des diagonal impliziten RK Verfahrens.

        Vergleiche `ExplicitRungeKutta.step`.
        """
        A, b, c = self.A, self.b, self.c
        s = self.n_stages
        dydt = np.empty((y0.size, s))   # werden auf Papier `k` genannt.

        for i in range(s):
            #def F(y_i):
            #    dydt_i = rhs(t + c[i]*dt, y_i)
            #    y_star = y0 + dt*(np.sum(A[i,:i]*dydt[:,:i], axis=-1) + A[i,i]*dydt_i)
            #    return y_i - y_star

            #initial_guess = y0
            #dydt[:,i] = rhs(t + c[i]*dt, fsolve(F, initial_guess))

            def F(x):
                x_star = rhs(t + c[i]*dt, y0 + dt*(np.sum(A[i,:i]*dydt[:,:i], axis=-1) + A[i,i]*x))
                return x - x_star

            initial_guess = y0#dydt[:,i]
            dydt[:,i] = fsolve(F,initial_guess)
        y1 = y0 + dt*np.sum(b*dydt, axis=-1)
        return y1


class ImplicitRungeKutta(RungeKutta):
    """Voll implizites Runge-Kutta Verfahren."""

    def step2(self, rhs, y0, t, dt):
        """Berechnet einen Schritt des impliziten RK Verfahrens.

        Vergleiche `ExplicitRungeKutta.step`.
        """
        A, b, c = self.A, self.b, self.c
        s = self.n_stages

        def F(dydt):
            dydt_star = np.empty_like(dydt)
            for i in range(s):
                y_star = y0 + dt*np.sum(A[i,:]*dydt[:,:], axis=-1)
                dydt_star[:,i] =  rhs(t + c[i]*dt, y_star)

            return dydt - dydt_star


        initial_guess = np.empty((y0.size, s))
        for i in range(s):
            initial_guess[:,i] = rhs(t, y0)

        dydt = nicer_fsolve(F, initial_guess)

        y1 = y0 + dt*np.sum(b*dydt, axis=-1)
        return y1

    def step(self, rhs, y0, t, dt):
        """Berechnet einen Schritt des diagonal impliziten RK Verfahrens.

        Vergleiche `ExplicitRungeKutta.step`.
        """
        A, b, c = self.A, self.b, self.c
        s = self.n_stages

        print("hello")
            # Mein Vorschlag:
        def F(x):
            x_star=np.empty_like(x)
            for i in range(s):
                x_star[:,i] = rhs(t + c[i]*dt, y0 + dt*np.sum(A[i,:]*x[:,:], axis=-1))# + A[i,i]*x))
            return x - x_star
        
        initial_guess = np.empty((y0.size, s))
        for i in range(s):
            initial_guess[:,i] = rhs(t, y0)
        dydt = nicer_fsolve(F,initial_guess)

        y1 = y0 + dt*np.sum(b*dydt, axis=-1)
        return y1


    


def nicer_fsolve(F, initial_guess):
    """Wrapper für `scipy.fsolve`.

    Dies ist nützlich, wenn man fsolve für 2D oder 3D Arrays
    brauchen will.
    """
    shape = initial_guess.shape

    initial_guess = initial_guess.reshape((-1,))
    result = fsolve(lambda x: F(x.reshape(shape)).reshape((-1)), initial_guess)

    return result.reshape(shape)


# Das Ordnung 5 Schema aus `ode45`, aka Fehlberg.
def exrk_o5():
    A = [[ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
         [ 0.25,  0.0,  0.0,  0.0,  0.0,  0.0],
         [ 3.0/32.0,  9.0/32.0,  0.0,  0.0,  0.0,  0.0],
         [ 1932.0/2197.0,  -7200.0/2197.0,  7296.0/2197.0,  0.0,  0.0,  0.0],
         [ 439.0/216.0,  -8.0,  3680.0/513.0,  -845.0/4104,  0.0,  0.0],
         [ -8.0/27.0,  2.0,  -3544.0/2565.0,  1859.0/4104.0,  -11.0/40.0,  0.0]]
    b = [16.0/135.0,  0.0,  6656.0/12825.0,  28561.0/56430.0,  -9.0/50.0,  2.0/55.0]

    return ExplicitRungeKutta(np.array(A), np.array(b))

# Runge-Kutta Butcher scheme for implicit midpoint rule
def dirk_im():
    A = np.array([[0.5]])
    b = np.array([1.0])

    return SemiImplicitRungeKutta(A, b)

# Diagonal implizites RK schema mit 3 Stufen und Ordnung 4.
def dirk_o4():
    alpha = 2.0*np.cos(np.pi/18)/np.sqrt(3)
    alpha2 = alpha*alpha

    A = np.array([[(1.0+alpha)/2.0,  0.0,  0.0,],
                  [-alpha/2.0,  (1.0+alpha)/2.0,  0.0],
                  [1.0+alpha,  -(1.0+2.0*alpha),  (1.0+alpha)/2.0]])
    b = np.array([1.0/(6.0*alpha2),  1.0-1.0/(3.0*alpha2),  1.0/(6.0*alpha2)])

    return SemiImplicitRungeKutta(A, b)


# Runge-Kutta Gauss Collocation of Order 4
def firk_o4():
    A = np.array([[0.25,  0.25 - np.sqrt(3)/6.0],
                  [0.25 + np.sqrt(3)/6.0,  0.25]])

    b = np.array([0.5, 0.5])

    return ImplicitRungeKutta(A, b)

# Runge-Kutta Gauss Collocation of Order 6
def firk_o6():
    A = np.array([[ 5.0/36.0,  2.0/9-np.sqrt(15)/15.0,  5.0/36.0-np.sqrt(15)/30.0],
                  [ 5.0/36.0+np.sqrt(15)/24.0,  2.0/9,  5.0/36.0-np.sqrt(15)/24.0],
                  [ 5.0/36.0+np.sqrt(15)/30.0,  2.0/9+np.sqrt(15)/15.0,  5.0/36.0 ]])

    b = np.array([5.0/18.0,  4.0/9.0,  5.0/18.0])

    return ImplicitRungeKutta(A, b)

if __name__ == "__main__":
    solver = dirk_o4()
    f = lambda t, y: -y
    y0 = np.array([1.0, 2.0])

    t_end, n_steps = 1.0, 1000

    t, y = solver(f, y0, t_end, n_steps)
    print(y[-1,:] - y0*np.exp(-t_end))
