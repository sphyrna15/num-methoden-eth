from numpy import zeros, tile
from scipy.optimize import fsolve


class RungeKutta(object):

    def __init__(self, A, b, c):
        """Build a new Runge-Kutta instance given a Butcher Scheme

        Input:
        A:  Butcher matrix A of shape (s,s)
        b:  Butcher vector b of shape (s,)
        c:  Butcher vector c of shape (s,)
        """
        self._A = A.copy()
        self._b = b.copy()
        self._c = c.copy()
        self._s = A.shape[0]


    def set_rhs(self, f):
        """Set the right-hand-side

        Input:
        f:  The right-hand-side f(t, y(t))
        """
        self._f = f


    def set_iv(self, T0, IV):
        """Set the initial values at time T0

        Input:
        T0:  The initial time T0
        IV:  The initial value y(T0)
        """
        self._T0 = T0
        self._IV = IV.copy()
        self._d = IV.shape[0]


    def integrate(self, Tend, steps):
        r"""
        Integrate ODE with Runge-Kutta Method

        Input:
        Tend:   Endzeit
        steps:  Anzahl Schritte

        Output:
        t, u:  Zeit und Loesung
        """
        u = zeros((self._d+1, steps+1))
        #################################################################################
        #                                                                               #
        # TODO: Implementieren Sie hier die Zeitevolution mittels Runge-Kutta Verfahren #
        #                                                                               #
        # Hinweis: Rufen Sie die Klassen-Methode self._step geeignet auf.               #
        #                                                                               #
        ##########################         DONE         #################################
        
        h = (Tend - self._T0) / float(steps)
        u[-1,0] = self._T0
        u[:-1,0] = self._IV
        
        for k in range(steps):
            u[:, k+1] = self._step(u[:,k], h )         
            
        return u[-1,:], u[:-1,:]


    def _step(self, u, dt):
        r"""
        Makes a single Runge-Kutta step of size dt, starting from current solution u(t).

        Input:
        u:     Current solution u(t) := [y(t), t]
        dt:    Timestep size

        Output:
        unew:  New solution u(t+dt) := [y(t+dt), t+dt]
        """
        d = self._d
        s = self._s
        unew = zeros((d+1,))
        ##################################################################################
        #                                                                                #
        # TODO: Implementieren Sie hier einen einzelnen Schritt eines                    #
        #       allgemeinen impliziten Runge-Kutta Verfahrens                            #
        #                                                                                #
        # Hinweis: Das Butcher Schema ist in den Klassen-Variablen self._A, self._b      #
        #          und self._c gespeichert.                                              #
        #                                                                                #
        #######################       DONE           #####################################
        
        ynew = zeros((d, ))
        
        K = self._compute_k(u, dt)
        
        y = zeros((d, ))
        for i in range(s):
            y += self._b[i] * K[:, i]
        
        unew[:-1] = u[:-1] + dt*y
        unew[-1] = u[-1] + dt
        
        return unew


    def _compute_k(self, u, dt):
        r"""
        Compute the k_i vectors for i = 1, ..., s.

        Input:
        u:     Current solution u(t) := [y(t), t]
        dt:    Timestep size

        Output:
        K:     Matrix of shape (d,s) whose columns are the k_i
        """
        pass


class ExplicitRungeKutta(RungeKutta):

    def _compute_k(self, u, dt):
        d = self._d
        s = self._s
        K = zeros((d, s))

        for i in range(s):
            targ = self._c[i]
            yarg = zeros((d,))
            for j in range(i):
                yarg += self._A[i,j] * K[:,j]

            K[:,i] = self._f(u[-1] + dt*targ, u[:-1] + dt*yarg)

        return K


class SemiImplicitRungeKutta(RungeKutta):

    def _compute_k(self, u, dt):
        d = self._d
        s = self._s
        K = zeros((s, d))

        # TODO zu implementieren.

        return K


class ImplicitRungeKutta(RungeKutta):

    def _compute_k(self, u, dt):
        d = self._d
        s = self._s

        # TODO zu implementieren.
        
        K = zeros((s,d))
        
        def F(k):
            temp = zeros((s*d, ))
            y = zeros((d, ))
            for i in range(s):
                for j in range(s):
                    
                    y += dt*self._A[i, j] * k[j*d : (j+1)*d]
                
                temp[i*d : (i+1)*d] = self._f(u[-1] + self._c[i]*dt, u[:-1] + y)
            
            return k - temp
        
        init = zeros((d, s))
        for i in range(s):
            init[:, i] = self._f(u[-1], u[:-1])
        K = fsolve(F, init)
        
        return K.reshape((s, d)).T
        
        
        # def F(k):
        #     tmp = zeros((s*d,))
        #     for i in range(s):
        #         targ = self._c[i]
        #         yarg = zeros((d,))
        #         for j in range(s):
        #             yarg += self._A[i,j] * k[j*d:(j+1)*d]

        #         tmp[i*d:(i+1)*d] = self._f(u[-1] + dt*targ, u[:-1] + dt*yarg)
        #     return tmp - k

        # K0 = tile(self._f(u[-1], u[:-1]), s)
        # K = fsolve(F, K0)
        # K = K.reshape((s,d)).T

        # return K


    