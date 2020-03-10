""" Beispiel Code aus der Uebungsstunde Woche3 """


""" Adaptive Quadratur """

def adaptint(f, a, b, tol):
    I_MPR = MPR(f, a, b)  # Mittelpunktregel ausführen
    I_SR = SR(f, a, b)    # Simpsonregel ausführen
    if abs(I_MPR - I_SR) < tol:   # falls unterschiel kleiner als tolerance (tol)
        return I_SR               # bessere approx zurück geben -> simpson

    else:  # Rekursion um die Genauigkeit des Integrals auf dem Intervall bis zur gewünschten tolerance erhöhen
        I_1 = adaptint(f, a, (a+b)/2, tol/2)  # falls nicht erzeugen wir eine neuen Mittelpunk
        I_2 = adaptint(f, a, (a+b)/2, tol/2)  # so können wir immer genauer werde -> kleinere Intervall

        I = I_1 + I_2
        return I       # noch bessere Implementation im Skript


# Adaptive Quadratur auch für zsm-gesetzte Quadraturformeln:

# Siehe Code in 1.4.2, im Skript S.23


""" Gauss-Legendre Quadraturformeln """

# Golub Welsch ALgorithmus im Skript 1.3.3 auf S.17
