# -*- coding: utf-8 -*-

import numpy as np

# TODO Implementieren Sie `trapezoid` und `simpson`
from quadrature import simpson, trapezoid

def trapezoid2d(f, a, b, Nx, c, d, Ny):
    """Approximiert das Doppelintegral mit der Trapezregel.

    Input:
           f : Funktion f(x, y) welche integiert werden soll.
        a, b : untere/obere Grenze des Integrals nach 'dx'.
          Nx : Anzahl Teilintervalle des Integrals nach 'dy'.
        c, d : untere/obere Grenze des Integrals nach 'dy'.
          Nx : Anzahl Teilintervalle des Integrals nach 'dy'.
    """
    y, dy = np.linspace(c, d, Ny+1, retstep=True)

    # Definiere,
    # F(y) = int_a^b f(x,y) dx
    F = lambda y: trapezoid(lambda x: f(x, y), a, b, Nx)

    # und integriere F(y) mit der Trapezregel.

    # TODO implementieren Sie.
    F = np.vectorize(F)
    I = trapezoid(F, c, d, Ny)
    return I

def simpson2d(f, a, b, Nx, c, d, Ny):
    """Approximiert das Doppelintegral mit der Trapezregel.

    Input:
           f : Funktion f(x, y) welche integiert werden soll.

        a, b : untere, obere Grenze des Integrals nach 'dx'.
          Nx : Anzahl Teilintervalle für das Integral über 'x'.

        c, d : untere, obere Grenze des Integrals nach 'dy'.
          Ny : Anzahl Teilintervalle für das Integral über 'y'.
    """

    # TODO wie trapezoid2d.
    
    #Definiere erst F(y)
    y, dy = np.linspace(c, d, Ny+1, retstep = True)
    F = lambda y: simpson(lambda x: f(x, y), a, b, Nx)
    F = np.vectorize(F)
    I = simpson(F, c, d, Ny)
    return I

def report_error(xp, yp, exact):
    """Berechnet und printet den Fehler im Punkt (xp, yp). """

    f = lambda x, y: 1.0 / np.sqrt((x - xp) ** 2 + (y - yp) ** 2)

    trapez_approx = trapezoid2d(f, -1.0, 1.0, 128, -1.0, 1.0, 128)
    trapez_error = trapez_approx - exact

    simpson_approx = simpson2d(f, -1.0, 1.0, 32, -1.0, 1.0, 32)
    simpson_error = simpson_approx - exact

    print('!- (xp,yp) = ( {:.1f}, {:.1f}) ---------------'.format(xp, yp))
    print('trapez (error)    : {:.8e} ({:+.8e}) '.format(trapez_approx, trapez_error))
    print('simpson (error)   : {:.8e} ({:+.8e}) '.format(simpson_approx, simpson_error))
    print()

if __name__ == "__main__":
    report_error(2.0, 2.0, exact = 1.449394876268660)
    report_error(10.0, 10.0, exact = 0.283080070385743)
    report_error(20.0, 20.0, exact = 0.141450870624223)
