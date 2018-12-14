 
import numpy as np

A = np.matrix([[1.,0,0,0],
               [0,2,0,0],
               [0,0,3,0],
               [0,0,0,4]])
eigval, eigvec = np.linalg.eig(A)

print(eigval)
print(eigvec)

# Create a "guess vector"
v = np.copy(eigvec[:,1])
v[0] = v[0] - 1e-4
v[1] = v[1] + 1e-5
v[2] = v[2] - 1e-2
v[3] = v[3] + 1e-5
v = v / np.linalg.norm(v,2)
print(np.linalg.norm(v))

print('Initial v \n', v )
print('Initial overlaps: \n', v.transpose() * eigvec)


x = np.linalg.solve(A,2*v)
print(x)
for i in range(100):
    v = np.linalg.solve(A, v)
    v = v / np.linalg.norm(v)
    print(v)



exit(0)

for i in range(100):
    M = A + 100*v.transpose() * np.eye(4)
    # M = A + v * v.transpose()*100
    v = M * v
    v = v / np.linalg.norm(v)
    overlaps = v.transpose() * eigvec
    print('Iteration ', i)
    print(' M*v:\n',v)
    print(" overlaps", overlaps)



# print('2 A*v:\n',A*v)


# print('A*v:\n',A*v)