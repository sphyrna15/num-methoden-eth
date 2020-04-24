import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import scipy.optimize

def omega2():
    g = 9.81
    l = 0.6

    return g/l

def pendulum_rhs(t, y):
    """Berechnet die rechte Seite der Pendelgleichung."""
    dydt = np.array([y[1], -omega2()*np.sin(y[0])])
    return dydt

def simulate_pendulum(initial_angle, t_end):
    """Berechnet die approximative Lösung der Pendelgleichung.

    Berechnet eine genaue Approximation der Pendelgleichung zur Zeit `t_end`,
    mit den Anfangsbedingungen:
        y0[0] = initial_angle
        y0[1] = 0.0

    Der ODE-Löser ist `solve_ivp`.

        initial_angle : Skalar. Auslenkung zur Zeit t = 0.
                t_end : Endzeit der Simulation.

    """
    # TODO: Implementieren Sie die Funktion.
    t, y = np.array([0.0, 1.0]), np.array([0.0, 0.0]) # dummy values.

    return t, y

def solve_for_final_angle(initial_angle, t_end):
    """Berechnet die Auslenkung des Pendels zur Zeit `t_end`."""
    # TODO: Implementieren Sie diese Funktion.
    return 0.

def simpson(f, a, b, N):
    """Berechnet das Integral von `f` über [a, b].

    Das Intgral wird mit Simpson Quadratur auf `N` Teilintervallen gelöst.

         f : Funktion welche integriert wird. Signatur: f(x), wobei x ein
             1D Array ist.
      a, b : Untere / obere Integrationsgrenze.
         N : Anzahl Teilintervalle.
    """
    # TODO: Berechnen Sie das Integral von f über [a, b].
    return 0.


def agm(x, y, n_iter=5):
    """Berechnet das arithmetisch-geometrische Mittel.

          x, y : zwei skalare.
        n_iter : Anzahl Iterationen des Verfahrens.
    """
    for n in range(n_iter):
        x, y = 0.5*(x + y), np.sqrt(x*y)

    return x

def plot_pendulum(t, y, c, phi):
    """Plotte Winkel und Winkelgeschwindigkeit als Funktion der Zeit.

          t : 1D Array. Zeit.
          y : 2D Array. Lösung der Pendelgleichung.
              Achse 0: Zeit.
              Achse 1: Winkel, Winkelgeschindigkeit
          c : Farbe (color) des Graphen, 'b', 'r', 'g', 'c', etc.
        phi : Anfangsauslenkung des Pendels.
    """

    phi_str = "{:1.2f}".format(phi/(0.5*np.pi))
    phi_label = r'$\phi \, | \, \phi_0 = ' + phi_str + r'\, \frac{\pi}{2}$'
    dphidt_label = r'$\frac{d\phi}{dt} \, | \, \phi_0 = ' + phi_str + r'\, \frac{\pi}{2}$'

    # TODO: Plotten Sie y[:,0] gegen t in der gefragten Farbe, mit label `phi_label`.
    # TODO: Plotten Sie y[:,1] gegen t in der gefragten Farbe, mit label `dphidt_label`.


def savefig(filename):
    """ Speichert den aktuellen matplotlib Plot als 'eps' und 'png'. """

    plt.savefig(filename + ".eps")
    plt.savefig(filename + ".png")

def ex2_c():
    t_end = 1.8

    phi_1, phi_2 = 0.0, 0.0  # TODO: fix
    t1, y1, = np.array([0.0, 0.0]), np.array([0.0, 0.0]) # TODO: fix
    t2, y2, = np.array([0.0, 0.0]), np.array([0.0, 0.0]) # TODO: fix

    # TODO: Lösen Sie Aufgabe 2 c)

    plot_pendulum(t1, y1, 'b', phi_1)
    plot_pendulum(t2, y2, 'r', phi_2)

    plt.title('Lösungen der Pendelgleichung')
    plt.xlabel(r'$t$')
    plt.xticks(np.linspace(0, t_end, 8+1))
    plt.legend(loc='upper left')
    plt.grid(True)

    savefig("ex2_c")

def ex2_d():
    print("Aufgabe 2 d) -- F(phi_0) = phi(0.45)")
    t_end = 1.8

     # TODO: Lösen Sie Aufgabe 2 d).
    phi_star = scipy.optimize.fsolve(F, initial_guess)
    phi_star = phi_star[0]
    print("phi_star = {:f}\n".format(phi_star))

    t3, y3 = simulate_pendulum(phi_star, t_end)
    plot_pendulum(t3, y3, 'g', phi_star)
    plt.legend(loc='upper left')
    savefig("ex2_d")

def ex2_e():
    print("Aufgabe 2 e) -- F(phi_0) = T(phi_0) - 0.45")

    # TODO: Lösen Sie Aufgabe 2 e).

    phi_star = np.array([0.0])  # TODO: fix
    print("phi_star = {:f}\n".format(phi_star[0]))

def ex2_f():
    print("Aufgabe 2 f) -- F(phi_0) = G(agm)")

    # TODO: Lösen Sie Aufgabe 2 f).
    phi_star = np.array([0.0])  # TODO: fix
    print("phi_star = {:f}\n".format(phi_star[0]))

if __name__ == '__main__':
    ex2_c()
    ex2_d()
    ex2_e()
    ex2_f()
    plt.show()

