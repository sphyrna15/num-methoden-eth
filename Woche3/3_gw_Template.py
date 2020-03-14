from numpy import *
from matplotlib.pyplot import *
from numpy.linalg import eigh


def trapezoidal(f, a, b, N):
    r"""
    Trapezoidal quadrature of function f from a to b with N subintervals.

    f:     Function f(x)
    a, b:  Bounds of the integration interval
    N:     Number of subintervals (2 evaluations per interval, internal nodes merged)
    """
    I = 0.0
    ############################################################
    #                                                          #
    # Implementieren Sie hier die zusammengesetzte Trapezregel #

    x, h = linspace(a, b, N+1, retstep = True)
    I = h * (sum(f(x[1:-1])) + 0.5 * (f(x[0]) + f(x[-1])))
    return I


def simpson(f, a, b, N):
    r"""
    Simpson quadrature of function f from a to b with N subintervals.

    f:     Function f(x)
    a, b:  Bounds of the integration interval
    N:     Number of subintervals (3 evaluations per interval, internal nodes merged)
    """
    I = 0.0
    #############################################################
    #                                                           #
    # Implementieren Sie hier die zusammengesetzte Simpsonregel #

    x, h = linspace(a, b, N+1, retstep = True)
    I = h/6 * (f(a) + 2*sum(f(x[1:-1])) + 4*sum(f((x[:-1] + x[1:]) * 0.5)) + f(b))
    return I


def gw_legendre(n):
    r"""
    Compute nodes and weights for Gauss-Legendre quadrature
    with the Golub-Welsch algorithm.

    n:  Number of node-weight pairs
    """
    x = zeros(n)
    w = zeros(n)
    i = arange(n-1)
    b = (i+1.0) / sqrt(4.0*(i+1)**2 - 1.0)
    J = diag(b, -1) + diag(b, 1)
    x, ev = eigh(J)
    w = 2 * ev[0,:]**2
    return x, w


def legendre(f, a, b, n):
    r"""
    Gauss-Legendre quadrature of function f from a to b with n nodes.

    f:     Function f(x)
    a, b:  Bounds of the integration interval
    n:     Number of quadrature evaluation points
    """
    I = 0.0
    ##############################################################
    #                                                            #
    # Implementieren Sie hier eine Gauss-Legendre Quadraturregel #
    nodes, weights = gw_legendre(n)
    
    transformed_nodes = 0.5 * ((1 - nodes)*a + (1 + nodes)*b)
    I = 0.5 * (b - a) * sum(weights * f(transformed_nodes))
    return I


def composite_legendre(f, a, b, N, n):
    r"""
    Composite Gauss-Legendre quadrature of function f from a to b
    with N subintervals, each using n nodes.

    f:     Function f(x)
    a, b:  Bounds of the integration interval
    N:     Number of subintervals
    n:     Number of quadrature evaluation points per interval
    """
    I = 0.0
    #################################################
    #                                               #
    # Implementieren Sie hier eine zusammengesetzte #
    # Gauss-Legendre Quadraturregel                 #
    
    x, h = np.linspace(a, b, N+1, retstep = True)
    
    F = np.vectorize(legendre)
    I = np.sum(F(f, x[:-1], x[1:], n)) 
    
    return I


def part_c():
    r"""All code for sub-task d)
    """
    # Plot
    t = linspace(0, 1, 200)
    figure(figsize=(12,8))
    plot(t, f1(t), label=r"$f_1(x) = \frac{1}{1 + 5x^2}$")
    plot(t, f2(t), label=r"$f_2(x) = \sqrt{x}$")
    grid(True)
    legend(loc="center right")
    # savefig("funktionen.png")

    # Exakte Werte der Integrale
    I1ex = arctan(sqrt(5.0)) / sqrt(5.0)
    I2ex = 2.0 / 3.0

    # Trapez Regel
    # Anzahl Teilintervalle
    Ntr = arange(1, 100)
    # Anzahl Auswertungspunkte insgesamt
    ntr = Ntr + 1

    I1tr = array([trapezoidal(f1, 0, 1, k) for k in Ntr])
    I2tr = array([trapezoidal(f2, 0, 1, k) for k in Ntr])
    E1tr = ones_like(I1tr)
    E2tr = ones_like(I2tr)
    #####################################
    #                                   #
    # Berechnen Sie den Quadraturfehler #
    
    E1tr = abs(I1ex - I1tr)
    E2tr = abs(I2ex - I2tr)

    # Simpson Regel
    # Anzahl Teilintervalle
    Nsi = arange(1, 100)
    # Anzahl Auswertungspunkte insgesamt
    nsi = 2*Nsi + 1

    I1si = array([simpson(f1, 0, 1, k) for k in Nsi])
    I2si = array([simpson(f2, 0, 1, k) for k in Nsi])
    E1si = ones_like(I1si)
    E2si = ones_like(I2si)
    #####################################
    #                                   #
    # Berechnen Sie den Quadraturfehler #
    
    E1si = abs(I1ex - I1si)
    E2si = abs(I2ex - I2si)

    # Gauss-Legendre
    # Anzahl Auswertungspunkte insgesamt
    ngl = arange(1, 100)

    I1gl = array([legendre(f1, 0, 1, k) for k in ngl])
    I2gl = array([legendre(f2, 0, 1, k) for k in ngl])
    E1gl = ones_like(I1gl)
    E2gl = ones_like(I2gl)
    #####################################
    #                                   #
    # Berechnen Sie den Quadraturfehler #
    
    E1gl = abs(I1ex - I1gl)
    E2gl= abs(I2ex - I2gl)

    # Composite Gauss-Legendre
    # Anzahl Teilintervalle
    Ncgl = arange(1, 50)
    # Anzahl Auswertungspunkte pro Teilintervall
    Kcgl = 4
    # Anzahl Auswertungspunkte insgesamt
    ncgl = Kcgl * Ncgl

    I1cgl = array([composite_legendre(f1, 0, 1, k, Kcgl) for k in Ncgl])
    I2cgl = array([composite_legendre(f2, 0, 1, k, Kcgl) for k in Ncgl])
    E1cgl = ones_like(I1cgl)
    E2cgl = ones_like(I2cgl)
    #####################################
    #                                   #
    # Berechnen Sie den Quadraturfehler #
    
    E1cgl = abs(I1ex - I1cgl)
    E2cgl = abs(I2ex - I2cgl)

    # Konvergenzplot fuer f1
    figure(figsize=(12,8))
    loglog(ntr, E1tr, "-o", label="Composite Trapez")
    loglog(nsi, E1si, "-o", label="Composite Simpson")
    loglog(ngl, E1gl, "-o", label="Gauss Legendre")
    loglog(ncgl, E1cgl, "-o", label="Composite Gauss Legendre")
    loglog(nsi, (1.0*nsi)**-2, "-k", label=r"$n^{-2}$")
    loglog(nsi, (1.0*nsi)**-4, "-.k", label=r"$n^{-4}$")
    grid(True)
    xlabel(r"Anzahl Auswertungen von $f$")
    ylabel(r"Absoluter Fehler")
    legend(loc="lower left")
    title(r"Quadrature von $f_1(x) = \frac{1}{1 + 5x^2}$")
    # savefig("konvergenz_f1.png")

    # Konvergenzplot fuer f2
    figure(figsize=(12,8))
    loglog(ntr, E2tr, "-o", label="Composite Trapez")
    loglog(nsi, E2si, "-o", label="Composite Simpson")
    loglog(ngl, E2gl, "-o", label="Gauss Legendre")
    loglog(ncgl, E2cgl, "-o", label="Composite Gauss Legendre")
    loglog(nsi, (1.0*nsi)**-2, "-k", label=r"$n^{-2}$")
    loglog(nsi, (1.0*nsi)**-4, "-.k", label=r"$n^{-4}$")
    grid(True)
    xlabel(r"Anzahl Auswertungen von $f$")
    ylabel(r"Absoluter Fehler")
    legend(loc="lower left")
    title(r"Quadrature von $f_2(x) = \sqrt{x}$")
    # savefig("konvergenz_f2.png")


def adaptquad(f, ML, rtol=1e-6, abstol=1e-10):
    r"""Adaptive quadrature using trapezoid and simpson rules

    f:       Function f(x)
    ML:      List of meshes, finest last
    rtol:    relative tolerance for termination
    abstol:  absolute tolerance for termination, necessary in case the exact
             integral value = 0, which renders a relative tolerance meaningless.
    """
    I = 0.0
    ############################################################
    #                                                          #
    # Implementieren Sie hier eine h-adaptive Quadratur        #
    #                                                          #
    # Hinweis: Passen Sie den Code 7.4.2 aus dem Skript so an, #
    #          dass er alle erzeugten Gitter zurueck gibt.     #
    M = ML[-1]
    h = diff(M)                     # compute lengths of mesh intervals
    mp = 0.5*(M[:-1] + M[1:])       # compute midpoint positions
    fx = f(M); fm = f(mp)           # evaluate function at positions and midpoints
    
    trp_loc = h*(fx[:-1] + 2 * fm+fx[1:]) / 4   # local trapezoid rule
    simp_loc= h*(fx[:-1] + 4 * fm+fx[1:]) / 6   # local simpson rule
    I = sum(simp_loc)                       # use simpson rule value as intermediate approximation for integral value
    est_loc = abs(simp_loc - trp_loc)   # difference of values obtained from 
                                        # local composite trapezoidal rule and local simpson rule is used as an estimate
                                        # for the local quadrature error.
    
    err_tot = sum(est_loc)              # estimate for global error (sum moduli of local error contributions)
    
    # if estimated total error not below relative or absolute threshold, refine mesh
    if err_tot > rtol * abs(I) and err_tot > abstol:
        
        refcells = nonzero(est_loc > 0.9*sum(est_loc)/size(est_loc))[0]
        new_M = sort(append(M,mp[refcells]))
        ML.append(new_M)
        
        I, ML = adaptquad(f , ML, rtol, abstol) # add midpoints of intervalls with large error contributions, recurse.
    
    return I, ML


def part_d():
    r"""All code for sub-task d)
    """
    # Adaptive Berechnung von f1
    M0 = linspace(0, 1, 5)
    I, ML = adaptquad(f1, [M0])

    figure(figsize=(12,8))
    plot(M0, zeros_like(M0), "db")
    for i, M in enumerate(ML):
        plot(M, i*ones_like(M), "|r")
    grid(True)
    xlabel(r"$x$")
    ylabel(r"Mesh refinement level")
    # savefig("mesh_f1.png")

    figure(figsize=(12,8))
    plot(list(map(len, ML)), "-o")
    grid(True)
    xlabel(r"Mesh refinement level")
    ylabel(r"Mesh size")
    # savefig("meshsize_f1.png")

    # Adaptive Berechnung von f2
    M0 = linspace(0, 1, 5)
    I, ML = adaptquad(f2, [M0])

    figure(figsize=(12,8))
    plot(M0, zeros_like(M0), "db")
    for i, M in enumerate(ML):
        plot(M, i*ones_like(M), "|r")
    grid(True)
    xlabel(r"$x$")
    ylabel(r"Mesh refinement level")
    # savefig("mesh_f2.png")

    figure(figsize=(12,8))
    plot(list(map(len, ML)), "-o")
    grid(True)
    xlabel(r"Mesh refinement level")
    ylabel(r"Mesh size")
    # savefig("meshsize_f2.png")




if __name__ == "__main__":
    # Funktionen f1 und f2 als globale Variablen definiert
    f1 = lambda x: 1.0 / (1.0 + 5.0*x**2)
    f2 = lambda x: sqrt(x)

    # Execute code for sub-tasks
    part_c()
    part_d()
