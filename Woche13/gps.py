import numpy as np
from scipy.constants import c
from numpy.linalg import lstsq, norm
import matplotlib.pyplot as plt
from matplotlib.mlab import find

## Flughoehe der GPS-Satelliten
rho = 26570*1e3
## Erdradius
re = 6370*1e3

def gn(x, F, J, tol=1e-6, maxit=100):
    """
    Gauss Newton algorithm.
    Hinweis: Damit Unteraufgabe c) funktioniert sollte mindestens `xp, ind`
    zurueckgegeben werden.

    Keyword Arguments:
    x     -- Startwert
    F     -- Residuenvektor
    J     -- Jacobi-Matrix
    tol   -- (default 1e-6)
    maxit -- (default 100)

    Returns:
    xp    --- Nullstelle von F
    ind   --- [True|False] Konvergenzindikator
    it    --- Iterationszaehler
    """
    #####################################################################
    # TODO: Implementiere hier Nonlinear-Least-Squares via Gauss-Newton #
    #####################################################################
    pass


def F(x, X):
    """
    Jacobi Matrix

    Keyword Arguments:
    x --  [ tr, xr, yr, zr]
    X --  [ ts, xs, ys, zs]
    """
    #################################################
    # TODO: Implementiere hier den Residuenvektor F #
    #################################################
    pass


def J(x, X):
    """
    Jacobi Matrix

    Keyword Arguments:
    x --  [ tr, xr, yr, zr]
    X --  [ ts, xs, ys, zs]
    """
    ###########################################################
    # TODO: Implementiere hier die Jacobi-Matrix der Funktion #
    ###########################################################
    pass


def cart_to_longlat(x):
    """
    Keyword Arguments:
    x -- position vector

    Returns:
    longtitude --- in degrees
    long. dir. --- direction O => east, or W => west
    latitude   --- in degrees
    lat. dir.  --- direction N => northern hemisphere, S => southern hemisphere
    """
    r = norm(x)
    theta = np.arcsin(x[2]/r)
    phi = np.arctan2(x[1],x[0])

    return np.rad2deg(phi), 'O' if x[1] > 0 else 'W', np.rad2deg(abs(theta)), 'N' if theta > 0 else 'S'

## Rotations matrizen
def rot_y(alpha):
    return np.array([[ np.cos(alpha), 0, np.sin(alpha) ],
                   [   0            , 1, 0             ],
                   [  -np.sin(alpha), 0, np.cos(alpha) ]])
def rot_z(alpha):
    return np.array([[ np.cos(alpha), np.sin(alpha), 0 ],
                   [  -np.sin(alpha), np.cos(alpha), 0 ],
                   [   0,             0,             1 ]])

## Code fuer Unteraufgabe c)
def random_sampling(satellite_positions, dtvar, plot_axes):
    """
    Fehleranalyse
    Plottet die Fehler in der Positionsbestimmung
    fuer 1000 zufaellig ausgewaehlte Untermengen der sichtbaren Satelliten mit
    Normal(0, dtvar)-verteilten Fehlern in der Zeitmessung.

    Keyword Arguments:
    satellite_positions ---
    dtvar               --- Standardabweichung der Messfehler
    plot_axes           --- matplotlib.axes object
    """
    np.random.seed(0)
    ax = plot_axes

    poe = re*np.array([1, 1, 1])/np.sqrt(3)
    # Finde die sichtbaren Satelliten
    I = np.dot(satellite_positions-poe, poe) > 0
    satellites_visible = satellite_positions[I,:]

    distance = norm( poe - satellites_visible, axis=1)
    travel_time = distance/c

    # add noise to measured travel time
    travel_time_measured = travel_time + np.random.normal(scale=dtvar, size=len(travel_time))
    # measurements
    X = np.hstack( (travel_time_measured.reshape(-1,1), satellites_visible) )

    x0 = np.array([0, 0, 0, 0])
    results = {i: [] for i in range(4,13)}

    for i in range(1000):
        idx = find(np.random.rand(len(X)) < 0.5)
        Xs = X[idx,:]
        if len(Xs) < 4:
            continue
        try:
            vals = gn(x0, lambda x: F(x,Xs), lambda x: J(x,Xs), tol=1e-6, maxit=100)
            x = vals[0]
            has_converged = vals[1]
            if has_converged:
                results[len(Xs)].append( (idx, norm(x[1:]-poe)) )
        except IndexError:
            print("Es ist ein Fehler aufgetreten.")
            print("Die Gauss-Newton Methode `gn` sollte mindestens folgende Werte zurueckgeben: (x, has_converged)")
            print("wobei: x die Loesung und has_converged [True|False] sein muss.")
            print("Passe gegebenenfalls deine Implementation an und versuch es nochmal.")
            return

    for key in results:
        results[key] = sorted( results[key], key= lambda x: x[-1])

    for key, value in results.items():
            ax.semilogy( len(value)*[key], [ v[1] for v in value], 'x')

    plt.grid(True)

    ax.set_xlabel('Anz. Satelliten')
    ax.set_ylabel('Fehler in der Position [m]')

if __name__ == '__main__':
    ## Verteile die Satelliten in ihren Umlaufbahnen
    N = 4 # Anz. Satelliten pro Bahn
    phi = 2*np.pi/N
    j = np.arange(0,N).reshape(-1,1)
    p0 = np.hstack([np.cos(j*phi), np.sin(j*phi), 0*j])

    satellite_positions = []
    for i in range(6):
        phi = np.pi/7
        pos_local = np.dot( rot_z(phi), p0.T).T
        pos_local = np.dot( rot_y(np.deg2rad(55)), pos_local.T).T
        pos_local = np.dot( rot_z(i*2*np.pi/6), pos_local.T).T
        satellite_positions.append(pos_local)
    satellite_positions = rho*np.vstack(satellite_positions)

    ## TODO: point on earth (poe) auswaehlen:
    #poe = np.array([1, 1, 1 ])/np.sqrt(3)*re
    #poe = np.array([1, 1, 0 ])/np.sqrt(2)*re
    poe = np.array([0, 0, 1 ])*re

    ## visible satellites
    I = np.dot(satellite_positions-poe, poe) > 0
    satellites_visible = satellite_positions[I,:]

    distance = norm(poe - satellites_visible, axis=1)
    travel_time = distance/c
    np.random.seed(0)
    dt = 1e-8 # Standardabweichung: 10 Nanosekunden
    ## Fuege normalverteilte Fehler hinzu
    travel_time_measured = travel_time + np.random.normal(scale=dt, size=len(travel_time))

    ## Messwerte
    X = np.hstack((travel_time_measured.reshape(-1,1), satellites_visible))

    ## Startwert
    x0 = np.array([0, 0, 0, 0])
    x, has_converged, it = gn(x0, lambda x: F(x,X), lambda x: J(x,X), tol=1e-8, maxit=100)

    ## Ausgabe
    if has_converged:
        print('Gauss-Newton liefert die Loesung, nach %d Iterationen:' % it)
        print('d = %.5e s' % x[0])
        print('x = %.3f km' % (x[1]/1000))
        print('y = %.3f km' % (x[2]/1000))
        print('z = %.3f km' % (x[3]/1000))
        long, OE, lat, SN = cart_to_longlat(x[1:])
        print('Wir befinden uns bei %.2f %s, %.2f %s' % (long, OE, lat, SN ))
    else:
        print('You are on the dark side of the moon.')


    ## d)
    ##########################################################################
    ## TODO: Falls du das Gauss-Newton Verfahren in a) und b) richtig        #
    ## implementiert hast, kannst du die folgenden Zeilen ausfuehren und die #
    ## Plots studieren.                                                      #
    ##########################################################################

    #plt.figure(figsize=(9,4))
    #ax = plt.subplot(141)
    #dt = 1e-10
    #random_sampling(satellite_positions, dt, ax)
    #plt.title(r'$\varepsilon_t = %.2e$' % dt)
    #plt.ylim(1e-2, 1e4)
    #plt.xlim(3,9)
    
    #ax = plt.subplot(142)
    #dt = 1e-8
    #random_sampling(satellite_positions, dt, ax)
    #plt.title(r'$\varepsilon_t = %.2e$' % dt)
    #plt.ylim(1e-2, 1e4)
    #plt.xlim(3,9)
    
    #ax = plt.subplot(143)
    #dt = 1e-7
    #random_sampling(satellite_positions, dt, ax)
    #plt.title(r'$\varepsilon_t = %.2e$' % dt)
    #plt.ylim(1e-2, 1e4)
    #plt.xlim(3,9)

    #ax = plt.subplot(144)
    #dt = 1e-6
    #random_sampling(satellite_positions, dt, ax)
    #plt.title(r'$\varepsilon_t = %.2e$' % dt)
    #plt.ylim(1e-2, 1e4)
    #plt.xlim(3,9)

    #plt.tight_layout()
    #plt.savefig('error-analysis.pdf')
    #plt.show()
