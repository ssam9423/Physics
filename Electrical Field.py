#!/usr/bin/env python
# coding: utf-8

# In[10]:


import numpy as np
import matplotlib.pyplot as plt


a = 1.0
b = 2.0
q = 1.0
e0 = 8.854187817 * 10**(-12)

def efieldx(X, Y):
    q1 = (-q * (X-a) / (4 * np.pi * e0 * (np.sqrt((X-a)**2 + (Y)**2))**3))
    q2 = (-q * (X+a) / (4 * np.pi * e0 * (np.sqrt((X+a)**2 + (Y)**2))**3))
    q3 = (-q * X) / (4 * np.pi * e0 * (np.sqrt((X)**2 + (Y+b)**2))**3)
    q4 = (3*q * X) / (4 * np.pi * e0 * (np.sqrt((X)**2 + (Y-b)**2))**3)
    return (q1 + q2 + q3 + q4)


def efieldy(X, Y):
    q1 = (-q * Y) / (4 * np.pi * e0 * (np.sqrt((X-a)**2 + (Y)**2))**3)
    q2 = (-q * Y) / (4 * np.pi * e0 * (np.sqrt((X+a)**2 + (Y)**2))**3)
    q3 = (-q * (Y+b) / (4 * np.pi * e0 * (np.sqrt((X)**2 + (Y+b)**2))**3))
    q4 = (3*q * (Y-b) / (4 * np.pi * e0 * (np.sqrt((X)**2 + (Y-b)**2))**3))
    return (q1 + q2 + q3 + q4)


# 4 by 4
L=4
nx, ny = 100, 100
x = np.linspace(-L, L, nx)
y = np.linspace(-L, L, ny)
X, Y = np.meshgrid(x, y)

# demo of streamplot for a vector field
fig = plt.figure(figsize = (8,8))
plt.title('PHYS 301 Problem Set 3 Question 3a')
plt.streamplot(x, y, efieldx(X, Y), efieldy(X, Y), linewidth=1, density=2, arrowstyle='->', arrowsize=1.5);


# 10 by 10
L2=10
nx2, ny2 = 100, 100
x2 = np.linspace(-L2, L2, nx2)
y2 = np.linspace(-L2, L2, ny2)
X2, Y2 = np.meshgrid(x2, y2)


# demo of streamplot for a vector field
fig = plt.figure(figsize = (8,8))
plt.title('PHYS 301 Problem Set 3 Question 3a')
plt.streamplot(x2, y2, efieldx(X2, Y2), efieldy(X2, Y2), linewidth=1, density=2, arrowstyle='->', arrowsize=1.5);




#############################################################



def exaprox(X, Y):
    r = np.sqrt(X**2 + Y**2)
    return (3*q*b*Y*X)/(np.pi*e0*r**5)

def eyaprox(X, Y):
    r = np.sqrt(X**2 + Y**2)
    return (q*b*(Y**2-2*X**2))/(np.pi*e0*r**5)

# 10 by 10
L3=10
nx3, ny3 = 100, 100
x3 = np.linspace(-L3, L3, nx3)
y3 = np.linspace(-L3, L3, ny3)
X3, Y3 = np.meshgrid(x3, y3)


# demo of streamplot for a vector field
fig = plt.figure(figsize = (8,8))
plt.title('PHYS 301 Problem Set 3 Question 3c')
plt.streamplot(x3, y3, exaprox(X3, Y3), eyaprox(X3, Y3), linewidth=1, density=2, arrowstyle='->', arrowsize=1.5);




# In[ ]:




