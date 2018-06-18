import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import cmath
import random


fig1 = plt.figure()
ax1 = Axes3D(fig1)


xvec = np.array([ -0.0050000000000000, -0.0040000000000000, -0.0030000000000000, -0.0020000000000000, -0.0010000000000000, 0.0000000000000000, 0.0010000000000000, 0.0020000000000000, 0.0030000000000000, 0.0040000000000000, 0.0050000000000000])
fvec = np.array([  1.0063887059347989, 1.0051073487971904, 1.0038278019772811, 1.0025500631308597, 1.0012741299173711, 0.9999999999999947, 0.9987276710456232, 0.9974571407248380, 0.9961884067119458, 0.9949214666849326, 0.9936563183254821 ])
gvec = np.array([  0.9999774880628309 + 0.0063660966688644j, 0.9999855923379807 + 0.0050928913225968j,  0.9999918956803858 + 0.0038196766513370j, 0.9999963980770802 + 0.0025464549863184j, 0.9999990995187997 + 0.0012732286587876j, 0.9999999999999956-0.0000000000000013j, 0.9999990995188051-0.0012732286587900j, 0.9999963980770820-0.0025464549863200j, 0.9999918956803843-0.0038196766513387j, 0.9999855923379809-0.0050928913225987j, 0.9999774880628325-0.0063660966688662j ])
ginf = gvec[6]
lambda_g2 = 0.9999991073691104-0.0012732286686835j
lambda_id = 1.0000000078503122+0.0000000000000001


print("First try: 1 site, no root \n")
G = lambda_g2/lambda_id
Gsq = G*G.conjugate()
Gabs = Gsq**0.5
Gsqrt = np.sqrt(G)
Gsqrt2 = Gsq**0.25
print("G  = " , G)
print("K2 = " , np.log(Gabs))
print("K2 = " , np.log(G**2))
print("K2 = " , np.log(Gsq))
print("K2 = " , np.log(np.abs(G)**2))
print("K2 = " , np.log(Gsqrt))
print("K2 = " , np.log(Gsqrt2))
print("K2 = " , 0.0000000075762938)

print("Second try: 2-site root on lambda_g2 only \n")
lambda_g2 = lambda_g2**0.5
G = lambda_g2/lambda_id
Gsq = G*G.conjugate()
Gabs = Gsq**0.5
Gsqrt = np.sqrt(G)
Gsqrt2 = Gsq**0.25
print("G  = " , G)
print("K2 = " , np.log(Gabs))
print("K2 = " , np.log(G**2))
print("K2 = " , np.log(Gsq))
print("K2 = " , np.log(np.abs(G)**2))
print("K2 = " , np.log(Gsqrt))
print("K2 = " , np.log(Gsqrt2))
print("K2 = " , 0.0000000075762938)


print("Third try: 2-site root on both \n")
lambda_g2 = lambda_g2**0.5
G = lambda_g2/lambda_id
Gsq = G*G.conjugate()
Gabs = Gsq**0.5
Gsqrt = np.sqrt(G)
Gsqrt2 = Gsq**0.25
print("G  = " , G)
print("K2 = " , np.log(Gabs))
print("K2 = " , np.log(G**2))
print("K2 = " , np.log(Gsq))
print("K2 = " , np.log(np.abs(G)**2))
print("K2 = " , np.log(Gsqrt))
print("K2 = " , np.log(Gsqrt2))
print("K2 = " , 0.0000000075762938)





print("Fourth try: identity transfer matrix \n")
lambda_g2 = 0.9999991073691104-0.0012732286686835j
G = lambda_g2
Gsq = G*G.conjugate()
Gabs = Gsq**0.5
Gsqrt = np.sqrt(G)
Gsqrt2 = Gsq**0.25
print("G  = " , G)
print("K2 = " , np.log(Gabs))
print("K2 = " , np.log(G**2))
print("K2 = " , np.log(Gsq))
print("K2 = " , np.log(np.abs(G)**2))
print("K2 = " , np.log(Gsqrt))
print("K2 = " , np.log(Gsqrt2))
print("K2 = " , 0.0000000075762938)


print("Fifth try: identity transfer matrix with upper root \n")
lambda_g2 = lambda_g2**0.5
G = lambda_g2
Gsq = G*G.conjugate()
Gabs = Gsq**0.5
Gsqrt = np.sqrt(G)
Gsqrt2 = Gsq**0.25
print("G  = " , G)
print("K2 = " , np.log(Gabs))
print("K2 = " , np.log(G**2))
print("K2 = " , np.log(Gsq))
print("K2 = " , np.log(np.abs(G)**2))
print("K2 = " , np.log(Gsqrt))
print("K2 = " , np.log(Gsqrt2))
print("K2 = " , 0.0000000075762938)





print("kappa2 = ", np.log(ginf**2))


gvec2 = np.zeros(shape=len(gvec),  dtype=complex)
for i in range(len(gvec)):
    # gvec2[i] = gvec[i].real
    gvec2[i] = gvec[i] * cmath.exp(-cmath.phase(gvec[i])*1.0j)


dgdx = np.diff(gvec2.real,1)/np.diff(xvec, 1)
dx  = 0.5 * (xvec[:-1] + xvec[1:])

d2gdx2 = np.diff(dgdx,1)/np.diff(dx, 1)
dx2 = 0.5 * (dx[:-1] + dx[1:])


N = len(dx)
slope = ( N * sum(np.multiply(dgdx,dx)) - sum(dgdx)*sum(dx)) / (N * sum(np.multiply(dx,dx)) - sum(dx)*sum(dx))
print("G'' = " , slope)

# ax1.scatter(xvec, gvec.real)
# ax1.scatter(xvec, gvec.real, gvec.imag)
# ax1.scatter(xvec, gvec2.real, gvec2.imag)
# plt.figure(2)
# plt.scatter(xvec, fvec)
# plt.figure(3)
# plt.scatter(dx, dgdx)
# plt.figure(4)
#
# plt.scatter(dx2, d2gdx2)
#
# plt.show()