"""

Seien die Matrizen A, C und die Vektoren b,d wie unten.


Berechnen Sie die Loesung $x_0$ des lineares Ausgleichsproblems $Cx=d$.
Schlagen Sie eine andere Loesung $x_g$ des Problems $Cx=d$ vor
Berechnen Sie die Loesung $x_a$ des lineares Ausgleichsproblems $Ax=b$.

Berechnen Sie die Loesung $x_m$ des lineares Ausgleichsproblems $Ax=b$ mit der
linearen Nebenbedingung $Cx=d$. 

Geben Sie fuer $x in $\{x_m, x_0, x_g, x_a\}$ aus:
+ $x$ 
+ die Residuen $r_A(x) = |b-Ax|_2$ und $r_C(x) = |d-Cx|_2$
+ die Normen $|x|_2$


import array_to_latex as a2l
a2l.to_ltx(C, frmt = '{:6.1f}', arraytype = 'array')
a2l.to_ltx(A, frmt = '{:6d}', arraytype = 'array')

A = 
\begin{array}
     5.0 &   -1.0 &   -1.0 &     6.0 &     4.0 &     0.0\\
   -3.0 &     1.0 &     4.0 &   -7.0 &   -2.0 &   -3.0\\
     1.0 &     3.0 &   -4.0 &     5.0 &     4.0 &     7.0\\
     0.0 &     4.0 &   -1.0 &     1.0 &     4.0 &     5.0\\
     4.0 &     2.0 &     3.0 &     1.0 &     6.0 &   -1.0\\
     3.0 &   -3.0 &   -5.0 &     8.0 &     0.0 &     2.0\\
     0.0 &   -1.0 &   -4.0 &     4.0 &   -1.0 &     3.0\\
   -5.0 &     4.0 &   -3.0 &   -2.0 &   -1.0 &     7.0\\
     3.0 &     4.0 &   -3.0 &     6.0 &     7.0 &     7.0
\end{array}


b =
\begin{array}
   -4.0 &     1.0 &   -2.0 &     3.0 &     3.0 &     0.0 &   -1.0 &     3.0 &     1.0
\end{array}


C =
\begin{array}
     1.0 &     3.0 &   -2.0 &     3.0 &     8.0 &     0.0\\
   -3.0 &     0.0 &     0.0 &     1.0 &     9.0 &     4.0\\
   -2.0 &     3.0 &   -2.0 &     4.0 &    17.0 &     4.0
\end{array}

d=
\begin{array}
     1.0 &     2.0 &   -3.0
\end{array}

"""


from numpy.linalg import matrix_rank, solve, svd, lstsq, norm
import numpy as np

A = np.array([ [5, -1, -1, 6, 4, 0],
               [-3, 1, 4, -7, -2, -3],
               [1, 3, -4, 5, 4, 7],
               [0, 4, -1, 1, 4, 5],
               [4, 2, 3, 1, 6, -1],
               [3, -3, -5, 8, 0, 2],
               [0, -1, -4, 4, -1, 3],
               [-5, 4, -3, -2, -1, 7],
               [3, 4, -3, 6, 7, 7]
               ])
m, n = A.shape

b = np.array([-4, 1, -2, 3, 3, 0, -1, 3, 1])
Ab = np.hstack((A,b.reshape(9,1) ))
print(matrix_rank(A ), matrix_rank(Ab) )

C = np.array([ [1, 3, -2, 3, 8, 0],
               [ -3, 0, 0, 1, 9, 4],
               [ -2, 3, -2, 4, 17, 4]
               ])

p,n = C.shape

d = np.array([1, 2, -3])
Cd = np.hstack((C,d.reshape(3,1) ))

rc = matrix_rank(C)
print(rc, matrix_rank(Cd) )


#####################################################################
# TODO: Finde x0, xg xa und xm #
#####################################################################


