# -*- coding: utf-8 -*-

import numpy as np
import scipy.optimize

def integrate(method, rhs, y0, T, N):
    """ Integriert die ODE vorwärts in Zeit.

        method : Numerische Methode um einen Zeitschritt zu machen.
           rhs : Rechte Seite der ODE.
            y0 : 1D Array, Anfangswerte der ODE.
             T : Endzeit.
             N : Anzahl Zeitschritte.
     """
    y = np.empty((N+1,) + y0.shape)

    t0, dt = 0, T/N # TODO: implement
    y[0,...] = y0
    for i in range(0, N):
        y[i+1,...] = method(rhs, y[i,...], t0 + i*dt, dt) # TODO: implement done

    t = np.arange(N+1)*dt
    return t, y

def explicit_euler_step(rhs, y0, t0, dt):
    """ Ein Schritt des expliziten Eulerverfahrens.

        rhs : Rechte Seite der ODE.
         y0 : Die aktuellen Werte der approximativen Lösung.
         t0 : Aktuelle Zeit.
         dt : Länge des Zeitschritts.
    """
    return y0 + dt*rhs(t0, y0)  # TODO: implement done

def explicit_euler(rhs, y0, T, N):
    return integrate(explicit_euler_step, rhs, y0, T, N)

def implicit_euler_step(rhs, y0, t0, dt):
    """ Impliziter Eulerschritt.

        Da `scipy.optimize.fsolve` nur 1D Arrays akzeptiert, müssen die
        Werte der Approximativen Lösung der ODE 1D Arrays sein.

        rhs : Rechte Seite der ODE.
         y0 : 1D Array. Die Werte zur aktuellen Zeit.
         t0 : Die aktuelle Zeit.
         dt : Länge des aktuellen Zeitschritts.
    """
     # TODO: Implementieren Sie das nicht-lineare Gleichungssystem. done
    F = lambda y: y - y0 - dt * rhs(t0+dt, y)

    initial_guess = y0 + dt * rhs(t0, y0)  # TODO: implementieren Sie einen besseren Startwert done
    return scipy.optimize.fsolve(F, initial_guess)

def implicit_euler(rhs, y0, T, N):
    return integrate(implicit_euler_step, rhs, y0, T, N)

def implicit_mid_point_step(rhs, y0, t0, dt):
    """ Implizite Mittelpunktsregel.

    Siehe `implicit_euler_step` für Dokumentation der Argumente.
    """
    F = lambda y: y - y0 - dt * rhs(0.5*(t0 + (dt+t0)), 0.5*(y0 + y))  # TODO: vgl. `implicit_euler_step` done
    
    initial_guess = y0 + dt * rhs(t0, y0)
    return scipy.optimize.fsolve(F, initial_guess)

def implicit_mid_point(rhs, y0, T, N):
    return integrate(implicit_mid_point_step, rhs, y0, T, N)

def velocity_verlet_step(rhs, xv0, t0, dt):
    """ Ein Zeitschritt des velocity-Verlet Verfahrens.

        rhs : Für das velocity-Verlet geeignete rechte Seite, d.h. die
              Beschleunigung des Systems. Die Beschleunigun hängt nur vom
              Ort ab.

        xv0 : Die Werte zur Zeit `t0`. Die erste hälfte des Arrays soll den Ort
              beinhalten, die zweite die Geschwindigkeit.

         t0 : Aktuelle Zeit (in der ODE).

         dt : Aktuelle Länge des Zeitschritts.
    """
    xv0 = xv0.reshape((2, -1))
    xv1 = np.empty_like(xv0)
    x0, x1 = xv0[0,:], xv1[0,:]
    v0, v1 = xv0[1,:], xv1[1,:]

    # TODO: Updaten Sie die Werte in `xv1` mittels velocity_verlet.
    x1[:] = x0 + dt*v0 + ((dt**2)*0.5) * rhs(t0, x0)
    v1[:] = v0+ 0.5* dt*(rhs(t0, x0) + rhs(t0+dt, x1))

    return xv1.reshape(-1)

def velocity_verlet(rhs, y0, T, N):
    return integrate(velocity_verlet_step, rhs, y0, T, N)
