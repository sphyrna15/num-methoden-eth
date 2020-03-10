# -*- coding: utf-8 -*-

import timeit
import numpy as np
import matplotlib.pyplot as plt
from ode_solvers import explicit_euler, implicit_euler, implicit_mid_point, velocity_verlet
from ode45 import ode45

def omega2():
    l = 0.6
    g = 9.81
    return g/l

def pendulum_rhs(t, y):
    """Rechte Seite der Pendulum ODE.

       t : Aktuelle Zeit.
       y : 1D Array. Approximative Lösung der Pendelgleichung.
           y[0] : Winkel.
           y[1] : Winkelgeschwindigkeit.
    """
    # TODO: implement
    return 0.0

def pendulum_verlet_rhs(t, alpha):
    """ Rechte Seitde der Pendelgleichung für velocity-Verlet.

        t : aktuelle Zeit.
        y : Aproximativer Winkel zur Zeit `t`.
    """
    # TODO: implement
    return 0.0

def convergence_rate(errors, resolutions):
    """ Berechnet die Konvergenzrate. """
    poly = np.polyfit(np.log(resolutions), np.log(errors), deg=1)
    return -poly[0]

def ex2_d():
    T = 0.2
    y0 = np.array([0.3, 0.0])

    all_methods = [explicit_euler, implicit_euler, implicit_mid_point, velocity_verlet]
    all_rhs = 3*[pendulum_rhs] + [pendulum_verlet_rhs]
    resolutions = 2**np.arange(4, 11)

    _, y_exact = ode45(pendulum_rhs, (0.0, T), y0, reltol=1e-12)

    for method, rhs in zip(all_methods, all_rhs):
        error = np.empty(resolutions.size)
        for k, N in enumerate(resolutions):
            # TODO: Berechen Sie die Lösung und den Fehler
            error[k] = 1.0

        rate = convergence_rate(error, resolutions)
        print(method)
        print("rate: " + str(rate) + "\n")

def potential_energy(y):
    """ Berechnet die potenzielle Energie von `y`.

        y : 2D-Array. Approximative Lösung der Pendelgleichung.
            Achse 0: Zeitschritte
            Achse 1: Winkel & Geschwindigkeit.
    """
    # TODO: Berechnen Sie die potenzielle Energie
    return 0.0

def kinetic_energy(y):
    """ Berechnet die kinetische Energie von `y`.

        y : 2D-Array. Approximative Lösung der Pendelgleichung.
            Achse 0: Zeitschritte
            Achse 1: Winkel & Geschwindigkeit.
    """
    # TODO: Berechnen Sie die kinetische Energie
    return 0.0

def create_plots(t, y, filename):
    E_pot = potential_energy(y)
    E_kin = kinetic_energy(y)
    E_tot = E_pot + E_kin

    plt.plot(t, E_pot, label="$E_{pot}$")
    plt.plot(t, E_kin, label="$E_{kin}$")
    plt.plot(t, E_tot, label="$E_{pot}$")
    plt.legend()

    plt.ylabel("Energy")
    plt.xlabel("Time")

    plt.savefig(filename + ".eps")
    plt.savefig(filename + ".png")
    plt.show()

def ex2_e():
    T, N = 4.0, 500
    y0 = np.array([1.4855, 0.0])

    all_methods = [explicit_euler, implicit_euler, implicit_mid_point, velocity_verlet]
    all_rhs = 3*[pendulum_rhs] + [pendulum_verlet_rhs]
    all_filenames = ["energy_ee", "energy_ie", "energy_im", "energy_vv"]

    for method, rhs, filename in zip(all_methods, all_rhs, all_filenames):
        t, y = method(rhs, y0, T, N)
        create_plots(t, y, filename)

def ex2_f():
    # TODO Benutzen Sie timeit um die Laufzeit zu messen.



if __name__ == '__main__':
    ex2_d()
    ex2_e()
    ex2_f()
