# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt

def lennard_acceleration(xy):
    """ Berechnet die Bescheunigung in einem Lennard-Jones Potential.

        xy : 1D array. Die (x,y)-Koordinaten des Teilchens.
    """
    # TODO: implementieren Sie den Gradient des Lennard-Jones
    #   Potentials.
    
    x, y = xy[0], xy[1]
    r2 = (x**2) + (y**2)
    
    # Mit Wolframalpha finden wir für den Gradient
    
    rdotdot = xy * 24 * (2 / (r2**7) - (r2**3) / (r2**7))
    return rdotdot

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

    # TODO: Implementieren Sie das Einschritt Störmer-Verlet Verfahren. done
    
    # Formeln aus Skript
    
    # Startwerte
    xy[0,:] = xy0
    vtemp = v0 + dt/2 * rhs(xy0)
    
    #Formel aus Skript
    for k in range(n_steps):
        xy[k+1,:] = xy[k,:] + dt*vtemp
        vtemp += dt*rhs(xy[k+1])


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
    
    # Auch hier nutzen wir wieder einfach die formeln aus dem Skript
    
    #Startwerte
    xy[0,:] = xy0
    xy[1,:] = xy0 + dt*v0 + (dt**2)/2 * rhs(xy0)
    
    #Formel
    for k in range(1, n_steps):
        xy[k+1,:] = -xy[k-1,:] + 2*xy[k,:] + (dt**2)/2 * rhs(xy[k])
    
    return t, xy

def simulate_trajectories(solver, label):
    """Führt die Simulation in Aufgabe 1 b) für `solver` durch.

        solver : Löser der verwendet werden soll.
         label : Ein String um die Löser zu unterscheiden.
     """

    # TODO: Berechnen Sie die Trajektorien der gefragten Teilchen
    #   und plotten Sie die Trajektorien.
     
    bgrid = np.linspace(0, 3, 20)
    t_end = 15
    dt = 0.02
    n_steps = int(t_end/dt)
    
    for b in bgrid:
        xy0 = np.array([-10, b])
        v0 = np.array([1, 0])
        t, xy = solver(lennard_acceleration, xy0, v0, n_steps, t_end)
        plt.plot(xy[:,0], xy[:,1])

    plt.grid(True)
    plt.xlabel('$x$')
    plt.ylabel('$y$')

    plt.savefig(label + ".eps")
    plt.savefig(label + ".png")

# if __name__ == '__main__':
#     simulate_trajectories(einschritt_stoermer_verlet, "einschritt_stoermer_verlet")
#     plt.show()

#     simulate_trajectories(zweischritt_stoermer_verlet, "zweischritt_stoermer_verlet")
#     plt.show()

