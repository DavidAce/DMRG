from sympy import *

t1, t2, t3, t4 = symbols('t1 t2 t3 t4', real=True)
c = symbols('c', complex=True)
M = Matrix([ [t1,0,0,0], [0,t2,c,0], [0,conjugate(c),t3,0], [0,0,0,t4] ])
res = M.eigenvects(simplify=True)
V = [v[2][0] for v in res]
V = Matrix.hstack(V[0],V[1],V[2], V[3])
Mtest = V.inv() * M * V
print(Mtest)
