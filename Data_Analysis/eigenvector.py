import numpy as np
import scipy.sparse.linalg.eigen.arpack as arp
import scipy.linalg as linalg
np.set_printoptions(linewidth=320)

M = np.loadtxt('file');
e0, v0 = arp.eigs(M, k=1, which='LR', return_eigenvectors=True)
v0 = v0.real
e0 = e0.real
v0[abs(v0)<1e-6] = 0
print('e0: \n', e0)
print('v0: \n', v0.reshape([8,8]))
