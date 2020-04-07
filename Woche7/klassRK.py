from numpy import size, linspace, zeros


def klassRK(f, y0, T, N):
    u"""
    Wendet auf eine System erster Ordnung der Form y'=f(t,y) das klassische
    Runge-Kutta-Verfahren an

    INPUT:
    f   - lambda function mit Argumenten (t,y)
    y0  - Vektor der Anfangswerte zum Zeitpunkt t0
    T   - T=[t0 tEnd] zu approximierendes Intervall
    N   - Anzahl der Iterationsschritte

    OUTPUT
    t   - Diskretisierung des Zeitintervalls
    y   - Matrix der approximierten Funktionswerte
    """
    # Anzahl der Gleichungen = Dimension der Loesung
    d = size(y0)
    # Zeitwerte und Schrittweite
    t, h = linspace(T[0], T[1], N + 1, retstep=True)

    y = zeros((d, size(t)))      # Loesungsmatrix
    y[:, 0] = y0                 # Anfangswert

    # Loesung berechnen
    for i in range(N):
        y[:, i + 1] = rkStep(f, y[:, i], t[i], h)

    return t, y


def rkStep(f, y0, t0, h):
    # Berechnet einen RK-Schritt der Schrittweite h
    # TODO: Implement RK-Schritt
    
    k1 = f(t0, y0)
    k2 = f(t0 + 0.5*h, y0 + h*(0.5*k1))
    k3 = f(t0 + 0.5*h, y0 + h*(0.5*k2))
    k4 = f(t0 + h, y0 + h*1.*k3)
    
    y1 = y0 + h*((1/6.)*k1 + (1/3.)*k2 + (1/3.)*k3 + (1/6.)*k4)
    
    return y1
