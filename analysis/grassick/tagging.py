import matplotlib.pyplot as plt
import numpy as np
import math
import random

##################################################################################

#                                    TAGGING A B AS B

##################################################################################

# Data from ATLAS NOTE,ATLAS-CONF-2015-057, 2015
x = [30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330,
   350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590]

b = [0.68, 0.778, 0.8, 0.809, 0.81, 0.81, 0.808, 0.8, 0.79, 0.783, 0.778, 0.765, 0.76,
     0.747, 0.735, 0.716, 0.704, 0.69, 0.67, 0.645, 0.64, 0.615, 0.6, 0.56, 0.546, 0.56,
     0.544, 0.523, 0.515]

# Truncated data on which a polynomial fit is applied
x2 = [30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330,
   350, 370, 390, 410, 430, 450, 470]

b2 = [0.68, 0.778, 0.8, 0.809, 0.81, 0.81, 0.808, 0.8, 0.79, 0.783, 0.778, 0.765, 0.76,
     0.747, 0.735, 0.716, 0.704, 0.69, 0.67, 0.645, 0.64, 0.615, 0.6]

z = np.polyfit(x2, b2, deg=10)
plt.xlabel("Jet Transverse Momentum [GeV/]", size=22)
plt.ylabel("Efficiency",size=22)
plt.title("B-tagging Performance",size=24)

def h(x):
    return z[0]*x**10 + z[1]*x**9 + z[2]*x**8 + z[3]*x**7 + z[4]*x**6 \
           + z[5]*x**5 + z[6]*x**4 + z[7]*x**3 + z[8]*x**2 + z[9]*x + z[10]



# Exponential tail
def f(x):
    return math.exp(-0.0012675*x + 0.0842395661)

# Piecewise defined function combining the above
def g(x):
    if x <= 470:
        return h(x)
    else:
        return f(x)



i = 0
a = []
c = []

while i < 1500:
    a.append(g(i))
    c.append(i)
    i += 1
plt.plot(c, a)
plt.scatter(x2,b2)
plt.show()




##################################################################################

#                                    MISTAGGING C AS B

##################################################################################

# Data from ATLAS NOTE,ATLAS-CONF-2015-057, 2015
x = [0, 30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330,
   350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590]

c = [0, 0.2, 0.23, 0.23, 0.225, 0.22, 0.22, 0.215, 0.21, 0.21, 0.2, 0.2, 0.195, 0.197,
     0.193, 0.19, 0.17, 0.167, 0.168, 0.16, 0.155, 0.13, 0.125, 0.131, 0.107, 0.107,
     0.12, 0.095, 0.09, 0.105]

# Polynomial fit
z = np.polyfit(x, c, deg=10)
plt.xlabel("Jet Transverse Momentum [GeV/]",size=22)
plt.ylabel("Efficiency",size=22)
plt.title("C-tagging Performance",size=24)

def h(x):
    return z[0]*x**10 + z[1]*x**9 + z[2]*x**8 + z[3]*x**7 + z[4]*x**6 \
           + z[5]*x**5 + z[6]*x**4 + z[7]*x**3 + z[8]*x**2 + z[9]*x + z[10]

#End tail of data used for exponential fit
x2 = [170, 190, 210, 230, 250, 270, 290, 310, 330,
   350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590]

c2 = [0.21, 0.21, 0.2, 0.2, 0.195, 0.197, 0.193, 0.19, 0.17, 0.167, 0.168, 0.16, 0.155,
      0.13, 0.125, 0.131, 0.107, 0.107, 0.12, 0.095, 0.09, 0.105]

j = 0
c3 = []

while j < len(c2):
    c3.append(math.log(c2[j]))
    j += 1

q = np.polyfit(x2, c3, deg=1)

def u(k):
    return math.exp(q[0]*k + q[1])

# Piecewise defined function combining the above
def f(x):
    if x < 435.334:
        return h(x)
    else:
        return u(x)

i = 0
a = []
d = []

while i < 1400:
    a.append(f(i))
    d.append(i)
    i += 1

c4 = []
for o in x:
    c4.append(h(o))

plt.plot(d, a)
plt.scatter(x,c)
plt.show()


##################################################################################

#                                    MISTAGGING C AS B

##################################################################################



# Data from ATLAS NOTE,ATLAS-CONF-2015-057, 2015
x = [30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330,
   350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590]

L = [0.0075, 0.0072, 0.007, 0.0068, 0.0065, 0.0068, 0.007, 0.0073, 0.0076, 0.0082, 0.0083, 0.0093, 0.0094,
     0.011, 0.0113, 0.0107, 0.012, 0.013, 0.0125, 0.01, 0.011, 0.0088, 0.0103, 0.0124, 0.011, 0.009,
     0.013, 0.0125, 0.0135]

# Polynomial fit
z = np.polyfit(x, L, deg=10)
plt.xlabel("Jet Transverse Momentum [GeV/]", size=22)
plt.ylabel("Efficiency",size=22)
plt.title("Light-flavour-tagging Performance",size=24)

def h(x):
    return z[0]*x**10 + z[1]*x**9 + z[2]*x**8 + z[3]*x**7 + z[4]*x**6 \
           + z[5]*x**5 + z[6]*x**4 + z[7]*x**3 + z[8]*x**2 + z[9]*x + z[10]

# Truncate function at 300 GeV due to lack of reliable data
def f(x):
    if x <= 300:
        return h(x)
    else:
        return h(300)

i = 0
a = []
c = []

while i < 590:
    a.append(f(i))
    c.append(i)
    i += 1

plt.plot(c,a)
# Truncated data
x2 = [30, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330]
L2 = [0.0075, 0.0072, 0.007, 0.0068, 0.0065, 0.0068, 0.007, 0.0073, 0.0076, 0.0082, 0.0083, 0.0093, 0.0094,
     0.011, 0.0113, 0.0107]
plt.scatter(x2,L2)
plt.show()
