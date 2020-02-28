# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

# ===================================================================================

# Unteraufgabe a)

def mittelpunkt(f,a,b,N):
    x, h = np.linspace(a, b, N+1, retstep = True)
    I = h * np.sum(f(0.5 * (x[:-1] + x[1:])))
    return I

def trapezoid(f, a, b, N):
    x, h = np.linspace(a, b, N+1, retstep = True)
    I = h * (np.sum(f(x[1:-1])) + 0.5 * (f(x[0]) + f(x[-1])))
    return I

def simpson(f,a ,b ,N):
    x, h = np.linspace(a, b, N+1, retstep = True)
    I = h/6 * (f(a) + 2*np.sum(f(x[1:-1])) + 4*np.sum(f((x[:-1] + x[1:]) * 0.5)) + f(b)) 
    return I
    
    
# ===================================================================================

# ===================================================================================

# Unteraufgabe b)
    

def quadrature_err(quad_rule, f, exact):
    n_chunks = 2 ** np.arange(3, 11)
    err = np.zeros((8, ))
    
    for i,N in enumerate(n_chunks):
        approximation = quad_rule(f,0 ,1 , N)
        err[i] = exact - approximation
    
    return err, n_chunks

# ===================================================================================

def plot_convergence(n_evals, errors, labels, title):
    """Plottet einen Konvergenzplot in einem neuen Fenster.

    Der Vergleich zwischen Simpson- und Trapezregel soll 'fair' sein. Die
    Simpsonregel wertet `f` fast 2x öfter aus. Plotten Sie deswegen den Fehler
    gegen die Anzahl Funktionsaufrufe.

    Trendlinien der Ordnung 1.5, 2 und 4 werden ebenfalls geplottet.

    Input:
        n_evals : list of array
                  Anzahl Funktionsaufrufe pro Regel pro Auflösung.

          errors : list of array
                  Fehler pro Quadraturregel pro Auflösung.

          labels : list. Text der in der Legende erscheint.
          title : Titel des Plots.
    """

    plt.figure(figsize=(12,8))

    for k, n_eval in enumerate(n_evals):
        plt.loglog(n_eval, errors[k], "-o", label=labels[k])

    # Trendlinien
    plt.loglog(n_evals[0], (1.0*n_evals[0])**-1.5, ":k", label=r"$n^{-1.5}$")
    plt.loglog(n_evals[0], (1.0*n_evals[0])**-2, "-k", label=r"$n^{-2}$")
    plt.loglog(n_evals[0], (1.0*n_evals[0])**-4, "-.k", label=r"$n^{-4}$")

    plt.grid(True)
    plt.xlabel("Anzahl Auswertungen von $f$")
    plt.ylabel("Absoluter Fehler")
    plt.legend(loc="lower left")
    plt.title(title)

# ===================================================================================

# Unteraufgabe c)
    
f1 = lambda x: 1 / (1 + 5 * x**2)
f2 = lambda x: np.sqrt(x)

# ===================================================================================

def convergence_experiment(f, exact, title):  #add filename variable to store images
    errors_mp, n_chunks = quadrature_err(mittelpunkt, f, exact)
    n_mp = n_chunks
    errors_tr, n_chunks = quadrature_err(trapezoid, f, exact)
    n_tr = n_chunks
    errors_si, n_chunks = quadrature_err(simpson, f, exact)
    n_si = n_chunks

    n_evals = [n_mp, n_tr, n_si]
    errors = [errors_mp, errors_tr, errors_si]
    labels = ["Mittelpunkt", "Trapez", "Simpson"]

    plot_convergence(n_evals, errors, labels, title)
#    plt.savefig(filename + ".png")
#    plt.savefig(filename + ".pdf")


def ex1_c():
    # Konvergenzplot für `f1`.
    I1ex = np.arctan(np.sqrt(5.0)) / np.sqrt(5.0)
    title = r"Quadraturfehler für $f_1(x) = \frac{1}{1 + 5x^2}$"
    convergence_experiment(f1, I1ex, title)  #filepath "img/convergence_f1"

    # Konvergenzplot für `f2`.
    I2ex = 2.0 / 3.0
    title = r"Quadraturfehler für $f_2(x) = \sqrt{x}$"
    convergence_experiment(f2, I2ex, title)  #filepath "img/convergence_f1"

if __name__ == "__main__":
    ex1_c()
    plt.show()
