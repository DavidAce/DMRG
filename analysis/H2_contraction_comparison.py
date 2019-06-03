import numpy as np
import matplotlib.pyplot as plt
x = np.arange(8,128,8)
a = 5
d = 2

left_right_with_mpos =\
    x * (2*x**2 + 2*d**4 + a**2 + a**4) \
   +d * (2*x**2 + 2*d**2 + a**2 + a**4) \
   +a * (2*x**2 + 2*d**2 + a**2)

current  = \
    x * (2*x**2 + d**2 + d**6 + a**2 + a**4) \
   +d * (2*x**2 + d**2 + d**4 + 2*a**2) \
   +a * (2*x**2 + 2*d**2 + a**2)




plt.figure()

plt.plot(x,left_right_with_mpos, label='LR w MPOS')
plt.plot(x,current, label='current')
plt.yscale('log')
plt.legend()
plt.show()