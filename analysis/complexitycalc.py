
import numpy as np

a = 2
c = 3
x = 8
i = 50
mat = c**2 * x**2 * a**2 + c**2 * x**2 * a**4  +  c* x**4 * a**4




tens = c*x**3*a**2 + 2*x**2*c*a**2*(c*a**2 + a) + x**2*a**2*(x*c**2 + c)   # "Fastest yet"

tens2 = x**2*c**2 *a**2 + x**2*c**2*a**2*(c*a**2+a) + x**2*a**2*c*(x*c**2 +c) +  x**2*a**2*(x*c**2 +c) # " theta - ha - hb - L - R "


tens3 = c**3*a**4 + x**2*c**2*a**4 + x**2*a**2*c*(x*a**4 + a**3 + a) + x**2*a**2*(x*c**2 + c) # " ha - hb - L - theta - R"


ten = x**2 * c**2 * a**2 + x**2*c*a**4 # store L ha hb first
op = x**2*a**2*c*(x*a**4 + a**3 + a) + x**2*a**2*(x*c**2 + c)


ten2 = 2*x**2*c**2*a**2 # store L ha  together with hb R first
op2  = x*a**2*c*(x*a**2 + a) + x**2*a**2 *(x*a**2*c**2 + c*a**2 + a)

ten3 = c**3*a**4 #1 gång
op3  = x**2*a**2*c*(c*a**4+ a**3 + a) + x**2*a**2*(x*c**2 + c)

print("tens                       : ", tens)
print("tens2                      : ", tens2)
print("tens3                      : ", tens3)

print("build mat                  : ", mat)
print("multiply mat               : ", x**4 * a**4)
print("build ten                  : ", ten)
print("multiply ten               : ", op)
print("build ten2                 : ", ten2)
print("multiply ten2              : ", op2)
print("build ten3                 : ", ten3)
print("multiply ten3              : ", op3)
print("50 iterations of mat       : ", mat + i * x**4 * a**4)
print("50 iterations of tens      : ", tens * i)
print("50 iterations of tens2     : ", tens2 * i)
print("50 iterations of tens3     : ", tens3 * i)
print("50 iterations of ten+op    : ", ten + i*op )

print("50 iterations of ten2+op2  : ", ten2 + i*op2 )
print("50 iterations of ten3+op3  : ", ten3 + i*op3 )
