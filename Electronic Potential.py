#!/usr/bin/env python
# coding: utf-8

# # Assignment 1: 
# In first year E&M you derived the electric potential of a charged spherical shell. If you forgot the answer, you can find it in Griffiths example 2.7. Plot this potential as a function of distance $r$ from the centre of the spherical shell of radius $R=10cm$ that carries a total charge $Q= 1\,\mu C$ for $0<r<50\,cm$. Make sure to label the axes of the graph and indicate appropriate units.
# 

# In[23]:


import numpy as np
import matplotlib.pyplot as plt


def V_o(r):                     # Electric potential outside of spherical shell
    return Q/(r*4*np.pi*e_0)

def V_i(r):                     # Electric potential inside of spherical shell
    n = len(r)                  # Potential within the shell is constant
    v = [Q/(R*4*np.pi*e_0)] * n
    return v
 
    
R   = 0.1                       # Radius of the spherical shell in meters
Q   = 0.000001                  # Total charge in Coulombs
e_0 = 8.85418782 * 10**(-12)    # Electric constant
dr  = 0.001
r_f = 0.5                       # Upper bound of distance r


r_i = np.arange(dr, R, dr)      # list of dist within shell
r_o = np.arange(R+dr, r_f, dr)  # list of dist outside shell


plt.figure()
plt.plot(r_o, V_o(r_o), c='g')
plt.plot(r_i, V_i(r_i), c='g')
plt.axvline(x=R, linestyle='--', c='grey', label='Radius of Spherical Shell')
plt.legend()
plt.grid()
plt.axis([0,r_f,0,100000])
plt.xlabel('Distance (m)')
plt.ylabel('Electric Potential (V)')


# # Assignment 2:
# 
# Consider a charged sphere of radius $R=10\,cm$ with charge density  
# 
# $$\rho(r) = 0.5e^{-r^2} {\rm C/m^3}$$  
# 
# Plot the total charge enclosed by a Gaussian (i.e. spherical) surface of radius r as a function of r for $0<r<20\,cm$. Again don't forget to label the axes of the graph.

# In[21]:


import functools as ft
import math


def C_den(r):                      # Charge density * r^2 (from spherical integration)
    pi4 = 4 * np.pi                # Integration of phi and theta
    r2  = np.square(r)
    return 0.5*pi4*r2/np.exp(r2)   # Value of C(r) at r


R   = 0.1                          # Radius of the spherical shell in meters
dr  = 0.001
r_f = 0.2                          # Upper bound of distance r

r_i = np.arange(dr, R+dr, dr)      # list of dist within shell
r_o = np.arange(R+dr, r_f, dr)     # list of dist outside shell


# Charge inside shell = C_den * dr at each point in list r_i
c_i = np.array([])                 # list of total charges inside sphere

for i in r_i:
    r_curr = np.arange(dr, i, dr)  # list of dist from 0 to current r
    c_curr = sum(C_den(r_curr)*dr) # integration of C(r)dr
    c_i = np.append(c_i, c_curr)
    
# Charge outside shell = constant
c_o = [sum(C_den(np.arange(dr, R, dr))*dr)] * len(r_o)


plt.figure()
plt.plot(r_o, c_o, c='g')
plt.plot(r_i, c_i, c='g')
plt.axvline(x=R, linestyle='--', c='grey', label='Radius of Spherical Shell')
plt.legend()
plt.grid()
plt.axis([0,r_f,0,0.0025])
plt.xlabel('Distance (m)')
plt.ylabel('Total Charge (C)')


# In[ ]:




