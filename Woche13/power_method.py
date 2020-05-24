# -*- encoding: utf-8 -*-
import numpy as np
import numpy.linalg

import scipy
import scipy.linalg

# Dies hilft mit der Formatierung von Matrizen. Die Anzahl
# Nachkommastellen ist 3 und kann hier ----------------------.
# geändert werden.                                           v
np.set_printoptions(linewidth=200, formatter={'float': '{: 0.3f}'.format})

ew = np.array([100.0, 10.0, 12.0, 0.04, 0.234, 3.92, 72.0, 42.0, 77.0, 32.0])
n = ew.size
Q, _ = np.linalg.qr(np.random.random((n, n)))
A = np.dot(Q.transpose(), np.dot(np.diag(ew), Q))

# TODO Finden sie den grössten Eigenwert vom A mittels Potenzmethode.

# TODO Finden sie den kleinsten Eigenwert vom A mittels inverser Potenzmethode.

# TODO Finden sie den Eigenwert am nächsten bei 42.
