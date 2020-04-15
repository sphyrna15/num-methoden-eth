# -*- coding: utf-8 -*-
import numpy as np
import scipy.optimize
import matplotlib.pylab as plt

# Wir benutzen die Löser aus Serie 2.
from ode_solvers import implicit_euler, explicit_euler

def rhs(t, y):
    """Die rechte Seite der ODE in Aufgabe 1."""
    # TODO zu implementieren.
    return np.array([ 0., 0.])

def solution(t):
    """Berechnet die exakte Lösung von Aufgabe 1.

    Input:
        t : ndarray
            1D array mit den t-Werten an denen die Lösung berechnet werden
            soll.
    """
    y = np.empty((t.shape[0], 2))

    # TODO zu implementieren.

    return y

def plot_comparison(all_t, all_y, all_labels, ylabel, filename):
    """Plottet einen Vergleich einer Komponente der Lösungen.

    Input:
             all_t : list
                     Zeit Gitter der jeweiligen Lösung.

             all_y : list
                     Die zu plottende Variable

        all_labels : list
                     Label in der Legende

            ylabel : string
                     Beschriftung der y-Achse

          filename : string
                     Der Plot soll hier gespeichert werden (ohne
                     Dateiendung).
    """
    all_styles = ['b-', 'r-', 'g-', 'k-']

    # TODO zu implementieren.

def ex1_b():
    ic = np.array([1.0, 0.0])

    t1, y1 = explicit_euler(rhs, ic, 2.0, 999)    # h > 0.002 --> unstable
    t2, y2 = explicit_euler(rhs, ic, 2.0, 2000)   # h < 0.002 --> stable
    t3, y3 = implicit_euler(rhs, ic, 2.0, 500)    # unconditionally stable
    y = solution(t1)

    all_t = [t1, t2, t3, t1]
    all_y = [y1[:,0], y2[:,0], y3[:,0], y[:,0]]
    all_labels = ["Explicit, unstable", "Explicit, stable", "Implicit", "exact"]

    plot_comparison(all_t, all_y, all_labels, r"$y_{0}", "img/stiff_u")

    all_y = [y1[:,1], y2[:,1], y3[:,1], y[:,1]]
    plot_comparison(all_t, all_y, all_labels, r"$y_{1}$", "img/stiff_v")

    plt.show()

if __name__ == "__main__":
    ex1_b()
