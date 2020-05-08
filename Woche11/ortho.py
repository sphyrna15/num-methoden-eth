import numpy as np
import matplotlib.pyplot as plt

def GR(a):
    """
    Gram-Schmidt Verfahren

    Keyword Arguments:
    a -- (m x n) Matrix

    Returns: q, r
    q -- Matrix Q
    r -- Matrix R
    """
    # TODO implement GS
    # DONE
    
    m, n = a.shape
    Q = np.zeros((m,n))
    R = np.zeros((n,n))
    
    for j in range(n):
        col_j = a[:,j]
        
        for i in range(j):
            R[i,j] = np.dot(Q[:,i], col_j)
            col_j = col_j - R[i,j] * Q[:,i]
        
        R[j,j] = np.linalg.norm(col_j)
        Q[:,j] = col_j / R[j,j]
            
    
    return Q, R


def GRmod(a):
    """
    Modifiziertes Gram-Schmidt Verfahren
    Keyword Arguments:
    a -- (m x n) Matrix

    Returns: q, r
    q -- Matrix Q
    r -- Matrix R
    """
    # TODO implement modified GS
    # DONE
    
    m, n = a.shape
    Q = np.zeros((m,n))
    R = np.zeros((n,n))
    V = a
    
    for j in range(n):
        R[j,j] = np.linalg.norm(V[:,j])
        Q[:,j] = V[:,j] / R[j,j]
        
        for k in range(j+1, n):
            R[k,j] = np.dot(Q[:,j], V[:,k])
            V[:,k] = V[:,k] - R[k,j] * Q[:,j]
    
    return Q, R


# Matrix Definition
n = 50
m = 50
Z = np.zeros((m, n))
for i in range(m):
    for j in range(n):
        Z[i, j] = 1 + min(i, j)

# numpy QR-implementation (als Vergleich)
q1, r1 = np.linalg.qr(Z)
# Guete des Gram-Schmidt Verfahren
q2, r2 = GR(Z)
# Guete des modifizierten Gram-Schmidt Verfahren
q3, r3 = GRmod(Z)
print("numpys qr liefert:         %.10e" % (np.dot(q1.T, q1) - np.eye(n)).max())
print("Gram-Schmidt liefert:      %.10e" % (np.dot(q2.T, q2) - np.eye(n)).max())
print("mod. Gram-Schmidt liefert: %.10e" % (np.dot(q3.T, q3) - np.eye(n)).max())

def plot_quality(q, label, filename):
    plt.figure()
    plt.title(r'$\log(|Q^T Q - I|)$, ' + label)
    im = plt.imshow(np.log10(abs(np.dot(q.T,q)-np.eye(n))+1e-16),
                    vmin=-16, vmax=1, interpolation='nearest')
    plt.colorbar()
    # plt.savefig(filename)

plot_quality(q1, "numpy.qr", "img/qr-errors-qr.png")
plot_quality(q2, "GS", "img/qr-errors-gs.png")
plot_quality(q3, "mod. GS", "img/qr-errors-gsmod.png")

plt.show()
