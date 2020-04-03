#!/usr/bin/env python
from numpy import *
from scipy.special import airy
from scipy.optimize import fsolve
from matplotlib.pyplot import *


# Anfangswerte
iv = array([ 0.35502805388781723926,
            -0.25881940379280679841])

# Exakte Loesung
exact = lambda t: row_stack(airy(t)[:2])

# Intervall
#TODO: Implement
T = 0


###################
#  Teilaufgabe 1  #
###################

# Airy rechte Seite
    #TODO: Implement
rhs = lambda t, y: 0


####################
#  Teilaufgabe 2   #
####################

def integrate_EE(y0, tstart, tend, N):
    r"""Integrate ODE with explicit Euler method

    Input: y0     ... initial condition
           tstart ... start t
           tend   ... end   t
           N      ... number of steps

    Output: t ... variable
            y ... solution
    """
    #TODO: Implement

    return tstart, y0


def integrate_IE(y0, tstart, tend, N):
    r"""Integrate ODE with implicit Euler method

    Input: y0     ... initial condition
           tstart ... start t
           tend   ... end   t
           N      ... number of steps

    Output: t ... variable
            y ... solution
    """
    #TODO: Implement

    return tstart, y0


def integrate_EM(y0, tstart, tend, N):
    r"""Integrate ODE with explicit midpoint method

    Input: y0     ... initial condition
           tstart ... start t
           tend   ... end   t
           N      ... number of steps

    Output: t ... variable
            y ... solution
    """
    #TODO: Implement

    return tstart, y0


def integrate_IM(y0, tstart, tend, N):
    r"""Integrate ODE with implicit midpoint method

    Input: y0     ... initial condition
           tstart ... start t
           tend   ... end   t
           N      ... number of steps

    Output: t ... variable
            y ... solution
    """
    #TODO: Implement

    return tstart, y0


N = 10**3
t, y = integrate_EE(iv, 0.0, T, N)
figure()
plot(t, exact(t)[0], "r", label = "Exakt")
plot(t, y[0,:], "b", label = "ExpEuler")
grid(True)
xlabel("t")
ylabel("Ai(t)")
legend(loc="best")
savefig("img/ai_expl_euler.png")

t, y = integrate_IE(iv, 0.0, T, N)
figure()
plot(t, exact(t)[0], "r", label = "Exakt")
plot(t, y[0,:], "b", label = "ImpEuler")
grid(True)
xlabel("t")
ylabel("Ai(t)")
legend(loc="best")
savefig("img/ai_impl_euler.png")

t, y = integrate_EM(iv, 0.0, T, N)
figure()
plot(t, exact(t)[0], "r", label = "Exakt")
plot(t, y[0,:], "b", label = "ExpMidpoint")
grid(True)
xlabel("t")
ylabel("Ai(t)")
legend(loc="best")
savefig("img/ai_expl_midpoint.png")

t, y = integrate_IM(iv, 0.0, T, N)
figure()
plot(t, exact(t)[0], "r", label = "Exakt")
plot(t, y[0,:], "b", label = "ImpMidpoint")
grid(True)
xlabel("t")
ylabel("Ai(t)")
legend(loc="best")
savefig("img/ai_impl_midpoint.png")


####################
#  Teilaufgabe 3   #
####################

def RK_step(B, t, y, h):
    r"""
    Makes a single Runge-Kutta step of size h.

    B: Butcher Schema
    t: Current time t_i
    y: Current solution y_i
    h: Timestep size
    """
    #TODO: Implement
    
    return t, y


def RK(B, y0, tstart, tend, N):
    r"""
    Integrate the equation with a Runge-Kutta rule.

    Input: B      ... Butcher Schema
           y0     ... initial condition
           tstart ... start t
           tend   ... end   t
           N      ... number of steps

    Output: t ... variable
            y ... solution
    """
    
    #TODO: Implement
    return tstart, y0

#TODO: Implement B
B = array([[0, 0],
           [0, 0]])

t, y = RK(B, iv, 0.0, T, N)
figure()
plot(t, exact(t)[0], "r", label = "Exact")
plot(t, y[0,:], "b", label = "RK_ExpEuler")
grid(True)
xlabel("t")
ylabel("Ai(t)")
legend(loc="best")
savefig("img/ai_RK_expl_Euler.png")

#TODO: Implement B
B = array([[0, 0],
           [0, 0]])

t, y = RK(B, iv, 0.0, T, N)
figure()
plot(t, exact(t)[0], "r", label = "Exact")
plot(t, y[0,:], "b", label = "RK_ImpEuler")
grid(True)
xlabel("t")
ylabel("Ai(t)")
legend(loc="best")
savefig("img/ai_RK_impl_Euler.png")

#TODO: Implement B
B = array([[ 0, 0, 0],
           [0,0, 0],
           [ 0, 0, 0]])

t, y = RK(B, iv, 0.0, T, N)
figure()
plot(t, exact(t)[0], "r", label = "Exact")
plot(t, y[0,:], "b", label = "RK_ExpMidpoint")
grid(True)
xlabel("t")
ylabel("Ai(t)")
legend(loc="best")
savefig("img/ai_RK_expl_midpoint.png")

#TODO: Implement B
B = array([[0,0],
           [ 0, 0]])

t, y = RK(B, iv, 0.0, T, N)
figure()
plot(t, exact(t)[0], "r", label = "Exact")
plot(t, y[0,:], "b", label = "RK_ImpMidpoint")
grid(True)
xlabel("t")
ylabel("Ai(t)")
legend(loc="best")
savefig("img/ai_RK_impl_midpoint.png")

####################
#  Teilaufgabe 4   #
####################

def RK_38(y0, tstart, tend, N):
    r"""
    Integrate the equation with the 3/8 Runge-Kutta rule.

    Input:
           y0     ... initial condition
           tstart ... start t
           tend   ... end   t
           N      ... number of steps

    Output: t ... variable
            y ... solution
    """

    #TODO: Implement

    return tstart, y0



t, y = RK_38(iv, 0.0, T, N)
figure()
plot(t, exact(t)[0], "r", label = "Exact")
plot(t, y[0,:], "b", label = "RK_38")
grid(True)
xlabel("t")
ylabel("Ai(t)")
legend(loc="best")
savefig("img/ai_RK_38.png")
