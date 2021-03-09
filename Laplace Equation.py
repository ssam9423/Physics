#!/usr/bin/env python
# coding: utf-8

# # Assignment
# 
# Generalize the procedure described above to two dimensions and solve in a compuational domain of size $N=50$
# 
# 1. the Laplace equation $$\Delta V(r)=0$$
# with the boundary conditions $V(x,0)=V(x,N)=\sin(\pi x/N)$ and $V(0,y)=V(N,y)=\sin(\pi y/N)$.
# 
# 
# 2. the Poisson equation $$\Delta V(r)=\delta^2 (r)=\delta(x)\delta(y)$$
# i.e. place a unit charge at the origin. Set $V(x,y)=0$ everywhere on the boundary.
# 
# Present your answer by visualizing the 2D array $V(x,y)$. One quick way to do this is to use the matplotlib function ${\tt pcolor(V)}$, which creates a colour plot on a grid with the colour representing the local value of V.  
# 

# In[6]:


import numpy as np
import matplotlib.pyplot as plt

# Setting Variables

N = 50         # Number of Points
i = 750        # Number of iterations
q = 1          # Unit charge

lap = np.zeros((N,N))
poi = np.zeros((N,N))


# Boundary Conditions
for ind in range(0, N - 1):
    lap[ind][0]   = np.sin(np.pi*ind / N)
    lap[ind][N-1] = np.sin(np.pi*ind / N)
    lap[0][ind]   = np.sin(np.pi*ind / N)
    lap[N-1][ind] = np.sin(np.pi*ind / N)
        
mid = int(N/2)
poi[mid][mid] = q

    
# Iterations
for it in range(0,i):                     # Goes through i iterations
    for indy in range(1, len(lap)-1):     # y indicies
        for indx in range(1, len(lap)-1): # x indicies
            lap[indx][indy]= (lap[indx-1][indy]+lap[indx+1][indy]) / 2
            poi[indx][indy]= (poi[indx-1][indy]+poi[indx+1][indy]) / 2
    for indx in range(1, len(lap)-1):     # x indicies
        for indy in range(1, len(lap)-1): # y indicies            
            lap[indx][indy]= (lap[indx][indy-1]+lap[indx][indy+1]) / 2
            poi[indx][indy]= (poi[indx][indy-1]+poi[indx][indy+1]) / 2

            
# Plot
plt.figure(0)
plt.title('Laplace Equation')
plt.pcolor(lap)

plt.figure(1)
plt.title('Poisson Equation')
plt.pcolor(poi)
            


# In[ ]:




