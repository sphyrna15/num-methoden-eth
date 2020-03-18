# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt

def lennard_acceleration(xy):
    """ Berechnet die Bescheunigung in einem Lennard-Jones Potential.

        xy : 1D array. Die (x,y)-Koordinaten des Teilchens.
    """
    # TODO: implementieren Sie den Gradient des Lennard-Jones
    #   Potentials.

    return np.zeros_like(xy)

def einschritt_stoermer_verlet(rhs, xy0, v0, n_steps, t_end):
    """ Einschritt Störmer-Verlet für ein Teilchen.

    Input: rhs     ... Rechte Seite der ODE
           xy0     ... Anfangskoordinaten xy0 = [x0,y0]
           v0      ... Anfangsgeschwindigkeiten v0 = [vx0,vy0]
           n_steps ... Anzahl Zeitschritte
           t_end   ... Endzeit

    Output: t ... Zeit
            xy ... Trajektorien xy = [x,y]
    """

    t, dt = np.linspace(0, t_end, n_steps+1, retstep=True)

    xy = np.zeros((n_steps+1,2))
    v = np.zeros((n_steps+1,2))

    # TODO: Implementieren Sie das Einschritt Störmer-Verlet Verfahren.

    return t, xy

def zweischritt_stoermer_verlet(rhs, xy0, v0, n_steps, t_end):
    """ Zweischritt Störmer-Verlet für ein Teilchen.

    Dieses Zweischrittverfahren kommt ohne extra Speicher für die
    Geschwindigkeit aus und halb so vielen Aufrufe der rechten Seite aus.

    Input: rhs     ... Rechte Seite der ODE.
           xy0     ... Anfangskoordinaten xy0 = [x0,y0]
           v0      ... Anfangsgeschwindigkeiten v0 = [vx0,vy0]
           n_steps ... Anzahl Zeitschritte
           t_end   ... Endzeit

    Output: t ... Zeit
            xy ... Trajektorien xy = [x,y]

    """

    t, dt = np.linspace(0.0, t_end, n_steps+1, retstep=True)
    xy = np.zeros((n_steps+1, 2))

    # TODO: Implementieren Sie das Zweischritt Störmer-Verlet Verfahren.

    return t, xy

def simulate_trajectories(solver, label):
    """Führt die Simulation in Aufgabe 1 b) für `solver` durch.

        solver : Löser der verwendet werden soll.
         label : Ein String um die Löser zu unterscheiden.
     """

    # TODO: Berechnen Sie die Trajektorien der gefragten Teilchen
    #   und plotten Sie die Trajektorien.

    plt.grid(True)
    plt.xlabel('$x$')
    plt.ylabel('$y$')

    plt.savefig(label + ".eps")
    plt.savefig(label + ".png")

if __name__ == '__main__':
    simulate_trajectories(einschritt_stoermer_verlet, "einschritt_stoermer_verlet")
    plt.show()

    simulate_trajectories(zweischritt_stoermer_verlet, "zweischritt_stoermer_verlet")
    plt.show()

