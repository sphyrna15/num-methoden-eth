""" Uebungsstunde Woche 2 Beispiel Code """

import numy as np
import matplotlib.pyplot as plt
import scipy 


""" Mittelpunktsregel """

# Intervall [a,b]

def MPR(f, a, b, N):
    X, h= np.linspace(a, b, N+1, retstep = True)
    I = h * np.sum(f(0.5 * X[:-1] + X[1:]))
    return I


""" Trapezregel """



""" Simpsonregel """



""" Exp Konvergenz """

N_max = 1000            #Anz Teilintervalle
N = np.arange(N_max + 1)

err = np.zeros(N_max)

for i in range(N_max):
    err[i] = abs(I_exact - QF(f, a, b, N[i]))

plt.loglog(N, err)
plt.loglog(N, N**(-2))

p = - np.polyfit(np.log(N), np.log(err), 1)[0]


""" Mehrdimensionale Quadratur """ 

#Intervalle [a,b] und [c,d] in 2-dim

F = lambda y: QF(lambda x: f(x,y), a, b, N)
I = QF(F, c, d, Ny)

#Reference Lösung aus SciPy
I_exact = scipy.integrate.nquad(f, np.array([[a, b], [c, d]]))[0]




