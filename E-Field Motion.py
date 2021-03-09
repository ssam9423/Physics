#!/usr/bin/env python
# coding: utf-8

# # PHYS 301 Problem Set 2
# 
# Computational Part
# 
# Sammy Song | 94237328

# # Motion of a Charge Through a Charged Ring
# 
# Write down the formula for the electric field along the axis of a charged ring of radius $R$ that carries a total charge $-Q$. Imagine that a charge $Q=1\,{\rm \mu C}$ of mass $m=100\,g$ is released from a point $x_0$ on the axis at time $t=0$. Compute the trajectory of the charge assuming $R=1m$ and that the ring carries the opposite charge $-Q$ up to $t=1000\,s$ for $x_0=1\,m$ and $x_0=10\,m$. 
# 
# In order to do this, we need to numerically integrate Newton's equation of motion ${\bf F}=m{\bf a}$. This is a second order ordinary differential equation. Without getting into any details about numerical methods for integrating ODEs, all of them are essentially based on a Taylor expansion in time for a small time time increment $\Delta t$. A suitable algorithm for this particular kind is the "velocity Verlet" method. Let us denote position, velocity and acceleration at time $t$ with $x(t)$, $v(t)$ and $a(t)$ then these quantities after a small time step $\Delta t$ can be obtained as follows:
# 
# 
# $$v(t+\Delta t/2)=v(t)+ a(t)\Delta t/2$$
# $$x(t+\Delta t)=x(t)+ v(t+\Delta t/2)\Delta t$$
# $$v(t+\Delta t)=v(t+\Delta t/2)+ a(t+\Delta t)\Delta t/2$$
# 
# These equations are in a form that can be directly implemented in a computer program. First define all relevant constants and a function that returns the force or acceleration as a function of position $x(t)$, and initialize the position to $x(0)=x_0$ and the velocity to $v(0)=0$. Then advance $x(t)$ and $y(t)$ according to the algorithm above and store the trajectory in an array. Your answer should consist of a plot of this trajectory (position versus time) for the two initial positions given above. Comment on the differences between these trajectories. Make sure you integrate long enough to see the characteristic pattern of the motion in each case.
# 
# Hints: By considering the symmetry of the problem, you should notice that the charge moves in one dimension only. You can try out different values of the timestep $\Delta t$ (what happens qualitatively), but a possible value for the parameters given above is $\Delta t=1s$. Also, ignore gravity. If you have time, experiment with the paramenters of the problem (mass, charge, size of ring, etc.)
# 
# 

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rcParams
#rcParams.update({'font.size':12}) #optional to change default font size

# define constants
m = 0.1
R = 1
Q = 10**(-6)

t_i = 0
t_f = 1000
dt = 1

x_01 = 1
x_02 = 10

v_0 = 0
a_0 = 0

e_0 = 8.8541878128 * 10**(-12)


# set up arrays for position and time
time = np.arange(t_i, t_f, dt)
index = np.arange(len(time) - 1)



# set initial position and velocity
pos1 = np.array(x_01)
pos2 = np.array(x_02)
vel1 = np.array(v_0)
vel2 = np.array(v_0)



# define a function that returns the force or acceleration
def acc(arr, ind):
    numerator = -Q*Q* arr[ind]
    denomenator = 4*np.pi*e_0*(R**2+arr[ind]**2)**(3/2)
    force = numerator / denomenator
    return force / m



# integrate equation of motion

# For initial position 1:
for i in index:
    if i == 0:
        num = -Q*Q* pos1
        den = 4*np.pi*e_0*(R**2+pos1**2)**(3/2)
        f = num / den
        a = f/m
        
        vel_half = vel1 + a*dt/2
        pos_full = pos1 + vel_half * dt
        pos1 = np.append(pos1, pos_full)
        
        vel_full = vel_half + acc(pos1, i+1) * dt/2
        vel1 = np.append(vel1, vel_full)
        
    else:
        a = acc(pos1, i)
    
        vel_half = vel1[i] + a*dt/2
        pos_full = pos1[i] + vel_half * dt
        pos1 = np.append(pos1, pos_full)
    
        vel_full = vel_half + acc(pos1, i+1) * dt/2
        vel1 = np.append(vel1, vel_full)
        
# For initial position 2:
for i in index:
    if i == 0:
        num = -Q*Q* pos2
        den = 4*np.pi*e_0*(R**2+pos2**2)**(3/2)
        f = num / den
        a = f/m
        
        vel_half = vel2 + a*dt/2
        pos_full = pos2 + vel_half * dt
        pos2 = np.append(pos2, pos_full)
        
        vel_full = vel_half + acc(pos2, i+1) * dt/2
        vel2 = np.append(vel2, vel_full)
        
    else:
        a = acc(pos2, i)
    
        vel_half = vel2[i] + a*dt/2
        pos_full = pos2[i] + vel_half * dt
        pos2 = np.append(pos2, pos_full)
    
        vel_full = vel_half + acc(pos2, i+1) * dt/2
        vel2 = np.append(vel2, vel_full)
            

    
# plot trajectories
plt.figure()
plt.plot(time, pos1, c='r', label=('x_0 = '+str(x_01)+' m'))
plt.plot(time, pos2, c='b', label=('x_0 = '+str(x_02)+' m'))
plt.legend()
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Distance (m)')



#####################    COMMENTS    #####################

# The behaviour of both charges is similar to that of a mass
# on a spring. This isn't that much of a surprise as the charges
# are free to move through the ring of charges, and there is no
# dampening force to slow the charges enough to slow them to a halt.

# The charge with the larger initial distance from the ring has 
# both a larger amplitude and a larger period than the charge with
# the smaller initial distance from the ring. This makes sense
# as the force acting on the charge initially is small, however
# as the charge gets closer to the ring, it speeds up a lot to
# the point where it shoots past the ring and travels 10 m in the 
# opposite direction until the force of the ring charge pulls 
# enough on the charge to a velocity of zero and repeat the cycle
# over and over again. 

