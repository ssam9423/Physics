#!/usr/bin/env python
# coding: utf-8

# In[25]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -

v_max = 25
dv = 0.1

r = 90/2
rho = 1.2
area = np.pi*r**2

vel = np.arange(dv, 25, dv)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def a_power(vel):
    return 3000*np.exp(vel-10) / (1+np.exp(vel-10))

def t_power(vel):
    return (rho*area*vel**3 / (2*1000))

def efficiency(vel):
    return (a_power(vel) / t_power(vel))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -

plt.figure()
plt.plot(vel, a_power(vel), c='b')
plt.grid()
plt.axis([0,v_max,0,3100])
plt.title('Actual Power')
plt.xlabel('Wind Speed (m/s)')
plt.ylabel('Actual Power (kW)')

plt.savefig('actualp.pdf')


plt.figure()
plt.plot(vel, t_power(vel), c='g')
plt.grid()
plt.axis([0,v_max,0,t_power(v_max)+10])
plt.title('Theoretical Power')
plt.xlabel('Wind Speed (m/s)')
plt.ylabel('Theoretical Power (kW)')

plt.savefig('theoreticalp.pdf')

plt.figure()
plt.plot(vel, efficiency(vel), c='r')
plt.axvline(x=12, linestyle='--', c='grey', label='12 m/s')
plt.legend()
plt.grid()
plt.axis([0,v_max,0,0.6])
plt.title('Efficiency')
plt.xlabel('Wind Speed (m/s)')
plt.ylabel('Efficiency')

plt.savefig('efficiency.pdf')


index12 = np.where(vel == 12)
print(efficiency(vel)[index12])


# In[ ]:




