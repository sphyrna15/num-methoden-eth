# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from ode_solvers import explicit_euler, implicit_euler
from ode_solvers import implicit_mid_point, velocity_verlet

def nbody_pdot(t, q, m, G):
    """ Berechnet die zeitliche Ableitung des Impluses im N-Körper Problem.

        t : Aktuelle Zeit, wird nicht gebraucht.

        q : 2D-Array, (n_bodies, 3). Ort jedes der N Körper. Entlang der ersten
            Achse sind die verschiedenen Körper, entlang der zweiten die drei
            örtlichen Dimensionen.

        m : 1D Array. Masse jedes Objekts.

        G : Gravitationskonstante.
    """
    dpdt = np.empty_like(q)

    # TODO Berechnene sie die zeitliche Ableitung von `p`
    # wir nutzen die Formel aus der Aufgabe
    
    for i in range(q.shape[0]):
        # Berechne Norm mit Python Broadcasting
        norm_dq = np.linalg.norm((q - q[i,:]), axis = 1) ** 3
        norm_dq = norm_dq.reshape((-1,1))
        # print("norm shape = %s" %str(norm_dq.shape))
        
        # Berechne Summen Term aus Formel
        sum_term = m * m[i] * (q - q[i,:]) / norm_dq
        # print("sumterm shape = %s" % str(sum_term.shape))
        
        # Nutze Python Broadcasting um über all den ganzen Array zu summieren
        dpdt[i,:] = G * (np.sum(sum_term[:i,:], axis = 0) + np.sum(sum_term[i+1:,:], axis = 0))
        
    return dpdt

def nbody_rhs(t, y, m, G, shape):
    """ Rechte Seite des N-Körperproblems.

            t : Aktuelle Zeit.
            y : 1D Array der approximativen Lösung der ODE.
            m : Masses der Objekte.
            G : Gravitationskonstante.
        shape : (2, n_bodies, 3). Die 'natürliche' Form von y. Entlang der ersten
                Achse sind Ort und Implus, entlang der zweiten die Körper und
                entlang der dritten die 3 örtlichen Dimensionen.
    """
    y = y.reshape(shape)
    q, p = y[0,...], y[1,...]
    m = m.reshape((-1,1))
    # print(q, p)

    dydt = np.empty_like(y)

    # TODO: implement done
    
    # Wir setzen nur die Formeln für qdot, pdot ein und schreiben
    # sie als Vektoren übereienander für die Ableitung von y
    dydt[0,...] = 1/m * p
    dydt[1,...] = nbody_pdot(t, q, m, G)

    return dydt.reshape(-1)

def nbody_verlet_rhs(t, q, m, G, shape):
    """ Zeitliche Ableitung der Geschwindigkeit des N-Körperproblems.

            t : Aktuelle Zeit.
            q : 1D Array des approximativen Ortes der Teilchen.
            m : Masses der Teilchen.
            G : Gravitationskonstante.
        shape : (2, n_bodies, 3). Die 'natürliche' Form von y. Entlang der ersten
                Achse sind Ort und Geschwindigkeit, entlang der zweiten die
                Teilchen und entlang der dritten die 3 örtlichen Dimensionen.
    """
    q = q.reshape(shape[1:])

    # TODO Implementieren Sie die rechte Seite für velocity-Verlet. done 
    
    # Wir wollen hier die Ableitung der Geschwindigkeit - Impuls ist ja genau
    # gegeben als Impuls = Masse * Geschwindigkeit - also teile
    # Impulsableitung durch m
    verlet_rhs = nbody_pdot(t, q, m, G) / m
    return verlet_rhs.reshape(-1)

def plot_orbits(y, filename):
    """ Erstellt einen Plot der Bahnen aller Teilchen.

               y : 4D-Array. Approximative Lösung des N-Körperproblems.
                   Achse 0: Zeitschritte
                   Achse 1: Ort und Implus/Geschwindigkeit.
                   Achse 2: Verschiedene Teilchen.
                   Achse 3: Raumdimensionen.

        filename : Name der Datei in der die Plot gespeichert werden.
    """

    for k in range(y.shape[2]):
        plt.plot(y[:,0,k,0], y[:,0,k,1])

    # plt.savefig(filename + ".png")
    # plt.savefig(filename + ".eps")
    plt.show()

def nbody_simulation(y0, m, T, N, G, figure_basename):
    n_bodies = m.size
    shape = y0.shape
    y0 = y0.reshape(-1)

    all_methods = [explicit_euler, implicit_euler, implicit_mid_point]
    rhs = lambda t, y: nbody_rhs(t, y, m, G, shape)
    all_labels = ["explicit_euler", "implicit_euler", "implicit_mid_point"]

    for method, label in zip(all_methods, all_labels):
        t, y = method(rhs, y0, T, N)
        plot_orbits(y.reshape((-1,) + shape), figure_basename + "_" + label)

    # Velocity-Verlet benötigt (q,v) als Startwerte.
    y0 = y0.reshape(shape)
    y0[1,:,:] /= m
    y0 = y0.reshape(-1)

    t, y = velocity_verlet(lambda t, y: nbody_verlet_rhs(t, y, m, G, shape), y0, T, N)
    plot_orbits(y.reshape((-1,) + shape), figure_basename + "_verlet")

def ex3_c():
    # TODO: Setzen Sie Anfangswerte done
    m = np.array([[500], [1]])
    G = 1 
    
    # Wir brauchen Startwerte - aus den vorherigen Teilaufgaben ist bekannt,
    # dass y die Form (2, n_bodies, 3) haben muss, also initialisieren wir:
    y0 =  np.zeros((2, 2, 3))
    
    # Fast alle Anfangswerte bleiben Null, ausser
    # Körper 1 startet bei Position (2, 0, 0):
    y0[0, 1, 0] = 2
    # Köper 2 startet mit Impuls (0, konstante, 0):
    y0[1, 1, 1] = np.sqrt((G * m[0] / 2))

    # TODO: Setzen Sie Integrationsparameter done
    T, N = 3, 5000
    

    nbody_simulation(y0, m, T, N, G, "two_body")

def ex3_d():
    n_bodies = 3
    m = np.array([1.0, 1.0, 1.0])
    G = 1.0

    shape = (2, n_bodies, 3)
    y0 = np.zeros(shape)

    y0[0,0,:] = np.array([0.97000436, -0.24308753, 0.0])
    y0[1,0,:] = np.array([0.46620368, 0.43236573, 0.0])

    y0[0,1,:] = -y0[0,0,:]
    y0[1,1,:] = y0[1,0,:]

    y0[0,2,:] = 0.0
    y0[1,2,:] = np.array([-0.93240737, -0.86473146, 0.0])

    T, N = 5.0, 10000

    nbody_simulation(y0, m, T, N, G, "eight")


def planet_ic(use_momentum=True):
    msun = 1.00000597682
    qsun = np.array([0,0,0])
    vsun = np.array([0,0,0])

    mj = 0.00095486104043
    qj = np.array([-3.5023653, -3.8169847, -1.5507963])
    vj = np.array([0.00565429, -0.00412490, -0.00190589])

    ms = 0.000285583733151
    qs = np.array([9.0755314, -3.0458353, -1.6483708])
    vs = np.array([0.00168318, 0.00483525, 0.00192462])

    mu = 0.0000437273164546
    qu = np.array([8.3101420, -16.2901086, -7.2521278])
    vu = np.array([0.00354178, 0.00137102, 0.00055029])

    mn = 0.0000517759138449
    qn = np.array([11.4707666, -25.7294829, -10.8169456])
    vn = np.array([0.00288930, 0.00114527, 0.00039677])

    mp = 7.692307692307693e-09
    qp = np.array([-15.5387357, -25.2225594, -3.1902382])
    vp = np.array([0.00276725, -0.00170702, -0.00136504])

    masses = np.array([msun, mj, ms, mu, mn, mp])

    n_bodies = 6
    y0 = np.empty((2, n_bodies, 3))
    for k, q0 in enumerate([qsun, qj, qs, qu, qn, qp]):
        y0[0,k,:] = q0

    for k, (m, v0) in enumerate(zip(masses, [vsun, vj, vs, vu, vn, vp])):
        if use_momentum:
            y0[1,k,:] = m*v0
        else:
            y0[1,k,:] = v0

    return masses, y0

def ex3_e():
    m, y0 = planet_ic(use_momentum=True)
    shape = y0.shape

    G = 2.95912208286e-4
    T, N = 20000, 2000

    nbody_simulation(y0, m, T, N, G, "solar_system")


if __name__ == '__main__':
    ex3_c()
    ex3_d()
    ex3_e()
