import numpy as np
import matplotlib.pyplot as plt

c = 3;
a = 2;

def LTRM(x):
    return x**3*c*a**2 + x**3*c**2*a**2 + x**2*a**3*(c**3*a**3 + c*a**3 + a**2 + 1)

def LTMR(x):
    return x**3*c*a**2 + x**2*c*a**2*(c*a**4 + a**3 + a) + x**2*a**2*(x*c**2 + c)

def LTmmR(x):
    return x**3*c*a**2 + 2*x**2*a**3*c*(c*a+1) + x**2*a**2*(x*c**2 + c)



xarr = np.arange(2,8)
plt.plot(xarr,LTRM(xarr), label='LTRM')
plt.plot(xarr,LTMR(xarr), label='LTMR')
plt.plot(xarr,LTmmR(xarr), label='LTmmR')
plt.legend()
plt.show()