import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from itertools import combinations
import scipy.optimize as sci

#====================================================================
# Constants

num_p = 400          # Number of particles
num_f = 1000         # Number of frames
equil = 400          # Number of frames until eqillibrium

dt = 0.00002         # Time step
size = 0.00001       # Size of particles
col = 5              # Factor that multiplies size

mass = 2.672E-26     # Mass of particles
k_b = 1.38064852E-23 # Boltzmann Constant

x_min, x_max = 0, 1
y_min, y_max = 0, 1

color_p = ['cyan', 'yellow', 'red']
dot_size = np.array([5])

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Figure

fig, ax = plt.subplots()
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)

#====================================================================
# Define Functions

def indices(comb_arr, orig_arr, coll1=True):
    col1_ind = []
    col2_ind = []
    
    if coll1:
        for i in comb_arr[:, 0]:
            index = np.where(orig_arr == i)
            col1_ind.append(index[0][0])
        return col1_ind
    else:
        for j in comb_arr[:, 1]:
            index = np.where(orig_arr == j)
            col2_ind.append(index[0][0])
        return col2_ind
    
    
def update_point(num):
    global x_pos, y_pos, x_vel, y_vel
    
    # Bounds: Any particles that are at the bounds bounce off elastically
    x_ind = np.where((x_pos < x_min) | (x_pos > x_max))
    y_ind = np.where((y_pos < y_min) | (y_pos > y_max))
    x_vel[x_ind] = -x_vel[x_ind]
    y_vel[y_ind] = -y_vel[y_ind]
    
    # Collisions: Finds particles that collide and bounce off elastically
    x_comb = np.asarray(list(combinations(x_pos, 2)))
    y_comb = np.asarray(list(combinations(y_pos, 2)))

    x_diff = x_comb[:, 0] - x_comb[:, 1]
    y_diff = y_comb[:, 0] - y_comb[:, 1]
    p_dist2 = (x_diff)**2+(y_diff)**2
    
    # Indices of Collision Pairs
    ind = np.where(p_dist2 < col * size)
    
    # Indicies of particles within collision pairs
    coll_ind1 = indices(x_comb[ind], x_pos)
    coll_ind2 = indices(x_comb[ind], x_pos, coll1=False)

    x_vel1, x_vel2 = x_vel[coll_ind1], x_vel[coll_ind2]
    y_vel1, y_vel2 = y_vel[coll_ind1], y_vel[coll_ind2]
    
    x_pos1, x_pos2 = x_pos[coll_ind1], x_pos[coll_ind2]
    y_pos1, y_pos2 = y_pos[coll_ind1], y_pos[coll_ind2]
    
    # Elastic Collision: Conserve both total energy and momentum
    # Update velocities
    factor = ((x_vel1 - x_vel2) * (x_pos1 - x_pos2) +
              (y_vel1 - y_vel2) * (y_pos1 - y_pos2)) / \
             ((x_pos1 - x_pos2)**2 + (y_pos1 - y_pos2)**2)
             
    x_vel[coll_ind1] = x_vel1 - factor * (x_pos1 - x_pos2)
    x_vel[coll_ind2] = x_vel2 + factor * (x_pos1 - x_pos2)
    y_vel[coll_ind1] = y_vel1 - factor * (y_pos1 - y_pos2)
    y_vel[coll_ind2] = y_vel2 + factor * (y_pos1 - y_pos2)
    
    # Update positions
    x_pos += dt * x_vel
    y_pos += dt * y_vel
    
    data = np.stack((x_pos, y_pos), axis=-1)
    image.set_offsets(data)
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Set Initial Conditions

x_pos = np.random.random(num_p)
y_pos = np.random.random(num_p)
x_vel = -500. * np.ones(num_p)

x_vel[np.where(x_pos <= 0.5)] = -x_vel[np.where(x_pos <= 0.5)]
y_vel = np.zeros(num_p)
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create Plot to be Animated

image = ax.scatter(x_pos, y_pos)
image.set_sizes(dot_size)

color_base = np.empty_like(x_pos, dtype=str)
color_base[np.where(x_pos <= 0.5)] = color_p[0]
color_base[x_pos >= 0.5] = color_p[1]
color_base[0] = color_p[2]
image.set_facecolor(color_base)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Animation

ani = animation.FuncAnimation(fig, update_point, frames = num_f, interval=10, 
                               repeat=False)

#animation.FuncAnimation.save(ani, filename="IdealGas.mp4", fps = 30)

#====================================================================
# Plotting the velocity distribution once at equilibrium

# Runs through collisions until at equilibrium
for i in np.zeros(equil):
    update_point(1)

#Velocity and Energy Distribuitions
v_distr = np.sqrt(x_vel**2 + y_vel**2)
E_distr = 1/2 * v_distr**2 * mass

def vel_true(vel, Temp):
    f_vel_analytical = ((mass * vel) / (k_b * Temp)) * \
        np.exp(-1/2 * mass * vel**2 / (k_b * Temp))
    return f_vel_analytical

def E_true(energy, Temp):
    f_energy_analytical = 1 / (k_b * Temp) * \
        np.exp(-energy/(k_b * Temp))
    return f_energy_analytical


# Theoretical distribution of Velocity and Energy
fitvel = np.linspace(v_distr.min(), v_distr.max(), 50)
fitenergy = np.linspace(E_distr.min(), E_distr.max(), 50)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Velocity Histogram
plt.figure(2)
plt.subplot(211)

n, bins, patch = plt.hist(v_distr, 50, normed=1, facecolor='g', 
                          histtype="bar")
fitparam, fitcov = sci.curve_fit(vel_true, fitvel, n, p0=200)

plt.plot(fitvel, vel_true(fitvel, *fitparam), '-r', label="Analytical Fit")
plt.title("Histogram of Velocity Distributions")
plt.legend()
plt.grid(True)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Energy Histogram
plt.subplot(212)

n, bins, patch = plt.hist(E_distr, 50, normed=1, facecolor='b', 
                          histtype="bar")
fitparam, fitcov = sci.curve_fit(E_true, fitenergy, n, p0=200)

plt.plot(fitenergy, E_true(fitenergy, *fitparam), '-r', label="Analytical Fit")
plt.title("Histogram of Energy Distributions")
plt.legend()
plt.grid(True)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Saves Histograms

#plt.tight_layout()
#plt.savefig("distributions.pdf")
