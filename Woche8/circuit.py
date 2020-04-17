from numpy import *
from matplotlib.pyplot import *
from scipy.optimize import fsolve
import time
from rk import *


# Implizite Mittelpunktsregel

def IM(u0, f, N, T, t0=0):
    """Implicit Midpoint Rule

    Input:
    u0:    Initial values u(t0) := [x(t0), x'(t0)]
    f:     Rechte-Seite f(t, u)
    N:     Anzahl Schritte
    T:     Endzeit
    t0:    Anfangszeit

    Output:
    t, u:  Zeit und Loesung
    """
    t, h = linspace(t0, T, N, retstep=True)
    u = zeros((u0.size, N))
    ################################################################
    #                                                              #
    # TODO: Implementieren Sie hier die implizite Mittelpunksregel #
    #                                                              #
    ###############          DONE           ########################
    
    u[:,0] = u0
    
    for i, time in enumerate(t[:-1]):
        
        F = lambda uj: u[:,i] + h*f(time+0.5*h, 0.5*u[:,i]+0.5*uj) - uj
        u[:, i+1] = fsolve(F, u[:,i] + h*f(time, u[:, i]))
    
    return t, u


# Butcher Schema für implizite Mittelpunktsregel
im_A = array([[0.5]])
im_b = array([1.0])
im_c = array([0.5])



# Runge-Kutta Gauss Collocation of Order 6
bs_A = array([[ 5.0/36.0,               2.0/9-sqrt(15)/15.0, 5.0/36.0-sqrt(15)/30.0],
              [ 5.0/36.0+sqrt(15)/24.0, 2.0/9,               5.0/36.0-sqrt(15)/24.0],
              [ 5.0/36.0+sqrt(15)/30.0, 2.0/9+sqrt(15)/15.0, 5.0/36.0              ]])

bs_b = array([5.0/18.0, 4.0/9.0, 5.0/18.0])

bs_c = array([0.5-sqrt(15)/10, 0.5, 0.5+sqrt(15)/10])


def run_im():
    # Alle benoetigten Parameter sind als globale Variablen verfuegbar
    Ioim = array([0.0, 0.0])
    #######################################################################
    #                                                                     #
    # TODO: Loesen Sie die Gleichung mit der impliziten Mittelpunktsregel #
    #       und obigen Anfangswerten Ioim = [y(t0), y'(t0)]               #
    #                                                                     #
    ################       DONE           #################################
    
    start = time.perf_counter()
    t, y = IM(Ioim, f, nsteps, T)
    end = time.perf_counter()
    imtime = end - start
    
    return "im", t, y, imtime


def run_rk():
    # Alle benoetigten Parameter sind als globale Variablen verfuegbar
    #################################################################
    #                                                               #
    # TODO: Achtung, Klassen-Template in rk.py umbenennen damit der #
    #       obige Import korrekt funktioniert                       #
    #################################################################
    Iork = array([0.0, 0.0])
    
    ########################################################################
    #                                                                      #
    # TODO: Loesen Sie die Gleichung mit der Runge-Kutta Gauss-Collocation #
    #       und obigen Anfangswerten Iork = [y(t0), y'(t0)]. Benutzen Sie  #
    #       die RungeKutta Klasse: erstellen Sie eine Instanz und rufen    #
    #       Sie die Member-Methoden korrekt auf.                           #
    #                                                                      #
    ###############          DONE          #################################
    RK = ImplicitRungeKutta(bs_A, bs_b, bs_c)
    RK.set_rhs(f)
    RK.set_iv(0.0, Iork)
    
    start = time.perf_counter()
    t, y = RK.integrate(T, nsteps)
    end = time.perf_counter()
    rktime = end - start
    
    print("RungeKutta Verfahren benötigt %s Sekunden" % (rktime))
    
    return "rk", t, y, rktime


def run_ode45():
    # Alle benoetigten Parameter sind als globale Variablen verfuegbar
    from ode45 import ode45
    Io45 = array([0.0, 0.0])
    ################################################################
    #                                                              #
    # TODO: Loesen Sie die Gleichung mit der dem ode45 Integrator  #
    #       und obigen Anfangswerten Io45 = [y(t0), y'(t0)]        #
    #                                                              #
    ###############        DONE          ###########################
    
    start = time.perf_counter()
    t, y = ode45(f, [0.0, T], Io45, stat='on', reltol = 1e-8, abstol = 1e-8, initialstep = 2e-5)
    end = time.perf_counter()
    ode45time = end - start
    
    print("Ode45 Verfahren benötigt %s Sekunden" % (ode45time))
    
    return "ode45", t, y.T, ode45time




if __name__ == '__main__':
    # Bauteile
    R = 100e3  # 100 KOhm
    C = 10e-9  # 10 nF
    L = 50e-3  # 50 mH

    # Signal
    Vo = 5  # 5 V
    w = 50  # 50 Hz
    Vin = lambda t: Vo * sin(2*pi*w*t)

    # Diodenmodell
    Is = 1e-9   # 1 nA
    Ut = 25e-3  # 25 mV
    n = 1.0
    Id = lambda t, yp: Is * (exp((Vin(t) - L*yp)/(n*Ut)) - 1.0)

    # Gleichung  f(t, u) = [u'(t), u''(t)]
    f = lambda t, y: array([0.,
                            0.])

    #############################################################
    #                                                           #
    # TODO: Implementieren Sie hier die rechte Seite f(t, y(t)) #
    #                                                           #
    ###############      DONE          ##########################
    
    f = lambda t, y: array([y[1], Id(t,y[1]) / (L*C) - y[1] / (R*C) - y[0] / (L*C)])

    # Number of periods
    np = 1
    # Number nodes per period
    nn = 8000
    # Number of steps in total
    #nsteps = np*nn + 2*nn/4 + 1
    nsteps = 12001
    # Endzeit
    T = np*0.02 + 2*0.005

    print(("Anzahl der Perioden: %d" % np))
    print(("Anzahl Zeitschritte: %d" % nsteps))
    print(("Endzeit: %f" % T))

    methodname = ""
    t = zeros((nsteps,))
    y = zeros((2,nsteps))
    ############################################################################
    #                                                                          #
    # TODO: Waehlen Sie hier die gewuenschte Methode durch aus/einkommentieren #
    #                                                                          #
    ############################################################################
    # methodname, t, y , methodtime = run_im()
    methodname, t, y, methodtime = run_rk()
    # methodname, t, y, methodtime = run_ode45()
    print(methodname + " Verfahren benötigt %s Sekunden" % (methodtime))

    ###################################################
    #                                                 #
    # TODO: Berechnen Sie hier die fehlenden Groessen #
    #                                                 #
    ###########    DONE    ############################

    Vout = y[1,:] * L
    IR = Vout / R
    IC = C*L * f(t, y)[1,:]
    ID = Id(t, y[1,:])
    IL = y[0,:]


    figure()
    plot(t, Vin(t), label=r'$V_{in}(t)$')
    plot(t, Vout, label=r'$V_{out}(t)$')
    xlabel(r'$t \; [s]$')
    ylabel(r'$V \; [V]$')
    grid(True)
    ylim(-Vo - 1, Vo + 1)
    legend(loc='lower right')
    # savefig('vosc_'+methodname+'.pdf')

    figure()
    plot(t, Vin(t), label=r'$V_{in}(t)$')
    plot(t, Vout, label=r'$V_{out}(t)$')
    xlabel(r'$t \; [s]$')
    ylabel(r'$V \; [V]$')
    grid(True)
    xlim(0.016, 0.021)
    ylim(-Vo - 1, Vo + 1)
    legend(loc='upper right')
    # savefig('vosc_zoom_'+methodname+'.pdf')

    figure()
    plot(t, squeeze(ID), label=r'$I_D(t)$')
    plot(t, squeeze(IR), label=r'$I_R(t)$')
    plot(t, squeeze(IC), label=r'$I_C(t)$')
    plot(t, squeeze(IL), label=r'$I_L(t)$')
    grid(True)
    xlabel(r'$t \; [s]$')
    ylabel(r'$I \; [A]$')
    legend(loc='upper right')
    # savefig('iosc_'+methodname+'.pdf')

    figure()
    #plot(t, squeeze(ID), label=r'$I_D(t)$')
    plot(t, squeeze(IR), label=r'$I_R(t)$')
    plot(t, squeeze(IC), label=r'$I_C(t)$')
    plot(t, squeeze(IL), label=r'$I_L(t)$')
    grid(True)
    xlim(0.016, 0.021)
    ylim(-0.0025, 0.0025)
    xlabel(r'$t \; [s]$')
    ylabel(r'$I \; [A]$')
    legend(loc='upper right')
    # savefig('iosc_zoom_'+methodname+'.pdf')

    show()