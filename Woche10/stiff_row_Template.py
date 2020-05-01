import numpy as np
from numpy.linalg import solve, norm

import matplotlib.pyplot as plt

###################
# Unteraufgabe a) #
###################

def integrate(f, Jf, y0, dt, n_steps, row_step):
    n_vars = y0.shape[0]
    t = np.arange(n_steps+1)*dt
    y = np.empty((n_steps+1, n_vars))

    y[0,:] = y0
    for i in range(n_steps):
        y[i+1,:] = row_step(f, Jf, y[i,:], dt)

    return t, y

def row_2_step(f, Jf, y0, dt):
    """Rosenbrock-Wanner Methode der Ordnung 2

    Input:
         f : Die rechte Seite der ODE f(x).

        Jf : Jacobi Matrix J(x) der Funktion, `shape == (n, n)`.

        y0 : ndarray.
             Aktueller Wert der approximativen Loesung der ODE.

        dt : Schrittweite

    Output:
        y1 : Zeitpropagierter Wert y(t+h).
    """
    ####################################################
    #                                                  #
    # TODO: Implementieren Sie die ROW-2 Methode hier. #
    #                                                  #
    a = 1.0 / (2.0 + np.sqrt(2.0))
    
    n = y0.shape[0]
    I = np.identity(n)
    J = Jf(y0)
    leftside = I - a*dt*J
    
    right_1 = f(y0)
    k1 = solve(leftside, right_1)
    
    right_2 = f(y0+0.5*dt*k1) - a*dt*np.dot(J,k1)
    k2 = solve(leftside, right_2)

    return y0 + dt*k2


def row_2(f, Jf, y0, dt, n_steps):
    return integrate(f, Jf, y0, dt, n_steps, row_2_step)

def row_3_step(f, Jf, y0, dt):
    """Rosenbrock-Wanner Methode der Ordnung 3

    Input:
         f : Die rechte Seite der ODE f(x).

        Jf : Jacobi Matrix J(x) der Funktion, `shape == (n, n)`.

        y0 : ndarray.
             Aktueller Wert der approximativen Loesung der ODE.

        dt : Schrittweite

    Output:
        y1 : Zeitpropagierter Wert y(t+h).
    """
    ####################################################
    #                                                  #
    # TODO: Implementieren Sie die ROW-3 Methode hier. #
    #                                                  #
    a = 1.0 / (2.0 + np.sqrt(2.0))
    d31 = - (4.0 + np.sqrt(2.0)) / (2.0 + np.sqrt(2.0))
    d32 = (6.0 + np.sqrt(2.0)) / (2.0 + np.sqrt(2.0))
    
    n = y0.shape[0]
    I = np.identity(n)
    J = Jf(y0)
    leftside = I - a*dt
    
    right_1 = f(y0)
    k1 = solve(leftside, right_1)
    
    right_2 = f(y0+0.5*dt*k1) - a*dt*np.dot(J,k1)
    k2 = solve(leftside, right_2)

    right_3 = f(y0+dt*k2) - d31*dt*np.dot(J,k1) - d32*dt*np.dot(J,k2)
    k3 = solve(leftside, right_3)
    
    
    return y0 + dt/6.0*(k1 + 4.0*k2 + k3)


def row_3(f, Jf, y0, dt, n_steps):
    return integrate(f, Jf, y0, dt, n_steps, row_3_step)

###################
# Unteraufgabe b) #
###################

class LogisticODE:
    def __init__(self, l):
        self.l = l

    def rhs(self, y):
        l = self.l
        return l*y*(1.0 - y)

    def jacobian(self, y):
        l = self.l
        Jf = l - 2.0*l*y[0]
        return Jf.reshape((1, 1))

    def exact(self, y0, t):
        l = self.l
        sol = (y0*np.exp(l*t)) / (1 - y0 + y0*np.exp(l*t))
        return sol.reshape((-1, 1))

def aufgabe_b():
    print(" Aufgabe b)")
    # Logistic ODE
    c = 0.01
    l = 25

    ode = LogisticODE(l)
    f = ode.rhs
    Jf = ode.jacobian

    y0 = np.array([c])

    T = 2.0
    N = 100
    h = T/float(N)

    ##################################################
    #                                                #
    # TODO: Loesen Sie die logistische Gleichung mit #
    #       den ROW Methoden und plotten Sie die     #
    #       Loesung und den Fehler.                  #
    #                                                #
    c = 0.01
    lamb = 25
    T = 2.0
    N = 100
    dt = T / float(N)
    
    ode = LogisticODE(l)
    f = ode.rhs
    J = ode.jacobian
    y0 = np.array([c])
    
    t2, y2 = row_2(f, J, y0, dt, N)
    t3, y3 = row_3(f, J, y0, dt, N)
    t = np.linspace(0.0, T, N+1)
    y = ode.exact(y0, t)
    
    plt.figure()
    plt.plot(t, y, "-go", label="Exakte Lösung $y(t)$, $\lambda = %s$" % lamb)
    plt.plot(t2, y2, "-r", label="ROW 2 $y(t)$, $\lambda = %s$" % lamb)
    plt.plot(t3, y3, "b", label="ROW 3 $y(t)$, $\lambda = %s$" % lamb)
    plt.title("ROW Methoden für logistische DLG")
    plt.xlabel(r"Zeit $t$")
    plt.ylabel(r"Lösung $y(t)$")
    plt.legend(loc='best')
    plt.grid()
    plt.show()
    
    plt.figure()
    plt.plot(t, np.abs(y2 - y), "-b", label="error ROW 2")
    plt.plot(t, np.abs(y3 - y), "-g", label="error ROW 3")
    plt.xlabel(r"Zeit $t$")
    plt.ylabel(r"Absolutbetrag Fehler")
    plt.title("Fehler für ROW Methoden")
    plt.legend(loc='best')
    plt.grid()
    plt.show()
    
    return None
    
    
    
    


###################
# Unteraufgabe c) #
###################

def aufgabe_c():
    print(" Aufgabe c)")
    # Logistic ODE
    c = 0.01
    ode = LogisticODE(l = 10.0)
    f = ode.rhs
    Jf = ode.jacobian

    y0 = np.array([c])
    T = 2.0

        ##################################################
        #                                                #
        # TODO: Loesen Sie die logistische Gleichung mit #
        #       den ROW Methoden.                        #
        #                                                #
    steps = 2**np.arange(4,13)
    e2 = np.empty(steps.size)
    e3 = np.empty_like(e2)

    for k, n in enumerate(steps):
        t, h = np.linspace(0.0, T, n+1, retstep=True)
        t2, y2 = row_2(f, Jf, y0, h, n)
        t3, y3 = row_3(f, Jf, y0, h, n)

        exact = ode.exact(y0, T)

        e2[k] = np.abs(y2[-1,0] - exact[0])
        e3[k] = np.abs(y3[-1,0] - exact[0])


    plt.figure()
    plt.loglog(steps, e2, "r", label="ROW 2")
    plt.loglog(steps, e3, "b", label="ROW 3")
    plt.loglog(steps, 1e-5*(T/steps)**2, "-k", label="$n^{-2}$")
    plt.loglog(steps, 1e-5*(T/steps)**3, "--k", label="$n^{-3}$")
    plt.xlabel(r"Anzahl Schritte $N$")
    plt.ylabel(r"Absoluter Fehler für $T = %.1f$" % T)
    plt.legend(loc="best")
    plt.grid()

    plt.show()


###################
# Unteraufgabe d) #
###################

def odeintadapt(Psilow, Psihigh, t_end, y0, fy0=None, dt0=None, dtmin=None, reltol=1e-2, abstol=1e-4):
    """Adaptive Integrator for stiff ODEs and Systems
    based on two suitable methods of different order.

    Input:
    Psilow:   Integrator of low order
    Psihigh:  Integrator of high order
    t_end:    Endtime
    y0:       Initial value, ndarray with shape == (n,).
    fy0:      The value f(y0): R^n used to estimate initial timestep size
    dt0:       Initial timestep (optional)
    dtmin:     Minimal timestep (optional)
    reltol:   Relative tolerance (optional)
    abstol:   Absolute Tolerance (optional)

    Output:
           t : Accepted time-steps
           y : Solution at the time-steps
    rejected : Array of rejected time-steps.
    error_estimates : Estimated error for each accepted time-step.
    """
    # Heuristic choice of initial timestep size
    if dt0 is None:
        dt0 = t_end / (100.0*(norm(fy0) + 0.1))
    if dtmin is None:
        dtmin = dt0 / 10000.0

    # Initial values
    yi = np.copy(y0)
    ti = 0.0
    dt = dt0

    # Use dynamically growing memory, since the algorithm
    # is adaptive.
    t = [ti]
    y = [yi]
    rejected = []
    error_estimates = [0.0]

    # while ti < t_end and dt > dtmin:
    #     #########################################################
    #     #                                                       #
    #     # TODO: Implementieren Sie hier die adaptive Strategie. #
    #     #                                                       #
    #     #########################################################

    n = yi.shape[0]
    t = np.array(t).reshape(-1)
    y = np.array(y).reshape(-1, n)
    rejected = np.array(rejected)
    error_estimates = np.array(error_estimates)
    return t, y, rejected, error_estimates


###################
# Unteraufgabe e) #
###################

def aufgabe_e():
    print(" Aufgabe e)")
    # Logistic ODE
    c = 0.01
    l = 50.0
    ode = LogisticODE(l)
    f = ode.rhs
    Jf = ode.jacobian
    y0 = np.array([c])
    T = 2.0

    # Discrete Evolution Operators
    Psilow = lambda dt, y: row_2_step(f, Jf, y, dt)
    Psihigh = lambda dt, y: row_3_step(f, Jf, y, dt)

    # Adaptive timestepping
    ta, ya, rej, ee = odeintadapt(Psilow, Psihigh, T, y0, f(y0))
    nsteps = ta.size

    print("Es werden %d Zeitschritte benoetigt" % nsteps)

    ########################################################################
    #                                                                      #
    # TODO: Testen Sie die adaptive Methode an der logistischen Gleichung. #
    #       Berechnen Sie die Anzahl der Zeitschritte. Plotten Sie die     #
    #       Loesung und den Fehler.                                        #
    #                                                                      #
    ########################################################################


###################
# Unteraufgabe f) #
###################

def aufgabe_f():
    print(" Aufgabe f)")
    # Test case
    tau = 4*np.pi/3.0
    R = np.array([[np.sin(tau), np.cos(tau)], [-np.cos(tau), np.sin(tau)]])
    D = np.array([[-101.0, 0.0], [0.0, -1.0]])
    A = np.dot(R.T, np.dot(D, R))
    print(A)

    f = lambda y: np.dot(A, y)
    Jf = lambda y: A
    y0 = np.array([[1.0],
                   [1.0]])
    T = 1.0

    ###################################################################
    #                                                                 #
    # TODO: Loesen Sie das Gleichungssystem mit der adaptiven Methode #
    #       und Plotten Sie die Loesung.                              #
    #                                                                 #
    ###################################################################


###################
# Unteraufgabe g) #
###################

def aufgabe_g():
    print(" Aufgabe g)")
    # Logistic ODE with y squared
    l = 500.0

    def f(y):
        return l*y**2*(1-y**2)

    def Jf(y):
        J = l*(2*y-4*y**3)
        return J.reshape((1,1))

    y0 = np.array([0.01])
    T = 0.5

    hmin = 1.0
    nrsteps = 0.0

    ##############################################################
    #                                                            #
    # TODO: Loesen Sie die steife Gleichung und plotten Sie      #
    #       die Loesung sowie die Groesse der Zeitschritte gegen #
    #       die Zeit.                                            #
    #                                                            #
    #       Wie viele Zeitschritte benoetigt dieses Verfahren?   #
    #       Was ist der kleinste Zeitschritt?                    #
    #       Wie viele Zeitschritte dieser Groesse wuerde ein     #
    #       nicht-adapives Verfahren benoetigen?                 #
    #                                                            #
    ##############################################################

    print("Minimal steps size: %f" % hmin)
    print("Number of adaptive steps: %i" % n_steps)
    print("Number of non-adaptive steps: %.2f" % (T/hmin))




if __name__ == "__main__":
    aufgabe_b()
    aufgabe_c()
    aufgabe_e()
    aufgabe_f()
    aufgabe_g()
