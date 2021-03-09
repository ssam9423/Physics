import numpy as np
import matplotlib.pyplot as plt

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# PHYSICS 216 201: HOMEWORK 4 PROBLEM 3
# Samantha Song | 94237328

# Problem #3: A 10 kg pumpkin is launched at 100 m/s from a height of 4 m
#             with a altitude of (0.942 * 360 = 339.12) and a azimuth of 
#             (0.942 * 70 = 65.94) at an initial velocity of 100 m/s
# Notes:      Assume that the atmosphere's gas is stationary and ignore the
#             effects of the Earth's rotation

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Set Up: Constants and Euler Functions

gravity = g = 9.8       # meters / second^2
height = z = 4          # meters
H = 8 * 1000            # meters | Constant
R = 6370 * 1000         # meters | Earth's Radius
delta_t = 0.0001        # seconds
mass = m = 10           # kilograms
initial_vel = 100       # meters / seconds
ang_vel = 7.29*10**-5   # radians / seconds

student_factor = 0.942
altitude = student_factor * 70   # degrees
azimuth = student_factor * 360   # degrees
latitude = 49           # degrees North

# Convert into Radians
alt = altitude / 180 * np.pi
azi = azimuth / 180 * np.pi
lat = latitude / 180 * np.pi


# Angular Velocity Components
ang_cos = ang_vel * np.cos(lat)
ang_sin = ang_vel * np.sin(lat)

# Time
time = [0]

# Actual Initial Positions, Velocities, and Acceleration
x_pos = [0]
y_pos = [0]
z_pos = [4]

x_vel = [initial_vel * np.cos(alt) * np.sin(azi)]
y_vel = [initial_vel * np.cos(alt) * np.cos(azi)]
z_vel = [initial_vel * np.sin(alt)]

x_acc = [0]
y_acc = [0]
z_acc = [-g - (0.043 * np.exp(-z_pos[-1]/H)) * z_vel[-1] * initial_vel / mass]

# Euler Functions | Essentially the same function
def euler_velocity(t, v_0, a):            # Returns Updated Velocity
    return a * t + v_0

def euler_position(t, p_0, v):            # Returns Updated Position
    return v * t + p_0

def corr_velocity(t, v_0, a):             # Returns Corrected Velocity
    return a * t + v_0

def corr_position(t, p_0, v):             # Returns Corrected Position
    return v * t + p_0

def pred_corr(euler, corrected):          # Returns Predictor Corrected Value
    return (euler + corrected)/2

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Numerially Euler Integrates Until Pumpkin Hits the Ground
def projectile(height = 4, initial_vel = 200, north = 'False', drag = 0.042): 
    if (north == 'North'):
        print('North')
        x_vel = [0]
        y_vel = [initial_vel]
        z_vel = [0]
    elif (north == 'East'):
        x_vel = [initial_vel]
        y_vel = [0]
        z_vel = [0]
        print('East')
    else:
        x_vel = [initial_vel * np.cos(alt) * np.sin(azi)]
        y_vel = [initial_vel * np.cos(alt) * np.cos(azi)]
        z_vel = [initial_vel * np.sin(alt)]
        
    # Time
    time = [0]
    # Actual Initial Positions, Velocities, and Acceleration
    x_pos = [0]
    y_pos = [0]
    z_pos = [height]
    
    x_acc = [0]
    y_acc = [0]
    z_acc = [-9.8 * (1 + z_pos[-1]/R) ** -2 - (drag * np.exp(-z_pos[-1]/H)) * 
             z_vel[-1] * initial_vel / mass]
    
    while (z_pos[-1] > 0):
        # Updates Time, Drag, Grabitational Force, and Total Velocity
        time.append(time[-1] + 1.5 * delta_t) # Current time
        b = drag * np.exp(-z_pos[-1]/H)       # Current b value
        g = 9.8 * (1 + z_pos[-1]/R) ** -2     # Current gravitational force
        total_vel = np.sqrt(x_vel[-1]**2 + y_vel[-1]**2 + z_vel[-1]**2)
        
        # Euler Acceleration, Velocity, and Position
        x_acc_eul = -b * x_vel[-1] * total_vel / mass
        x_acc_eul -= 2 * z_vel[-1] * ang_cos 
        x_acc_eul += 2 * y_vel[-1] * ang_sin - x_pos[-1] * ang_vel**2
        y_acc_eul = -b * y_vel[-1] * total_vel / mass
        y_acc_eul -= 2 * x_vel[-1] * ang_sin
        y_acc_eul -= ang_sin * (z_pos[-1] * ang_cos - y_pos[-1] * ang_sin)
        z_acc_eul = -g - b * z_vel[-1] * total_vel / mass 
        z_acc_eul += 2 * x_vel[-1] * ang_cos
        z_acc_eul += ang_cos * (z_pos[-1] * ang_cos - y_pos[-1] * ang_sin)
        
        x_vel_eul = euler_velocity(delta_t, x_vel[-1], x_acc[-1])
        y_vel_eul = euler_velocity(delta_t, y_vel[-1], y_acc[-1])
        z_vel_eul = euler_velocity(delta_t, z_vel[-1], z_acc[-1])
        
        x_pos_eul = euler_position(delta_t, x_pos[-1], x_vel[-1])
        y_pos_eul = euler_position(delta_t, y_pos[-1], y_vel[-1])
        z_pos_eul = euler_position(delta_t, z_pos[-1], z_vel[-1])
        
        euler_vel = np.sqrt(x_vel_eul**2 + y_vel_eul**2 + z_vel_eul**2)
        
        # Corrected Acceleration, Velocity, and Position
        x_acc_cor = -b * x_vel_eul * euler_vel / mass
        x_acc_eul -= 2 * z_vel_eul * ang_cos 
        x_acc_eul += 2 * y_vel_eul * ang_sin - x_pos_eul * ang_vel**2
        y_acc_cor = -b * y_vel_eul * euler_vel / mass
        y_acc_eul -= 2 * x_vel_eul * ang_sin
        y_acc_eul -= ang_sin * (z_pos_eul * ang_cos - y_pos_eul * ang_sin)
        z_acc_cor = -g - b * z_vel_eul * euler_vel / mass
        z_acc_cor += 2 * x_vel_eul * ang_cos
        z_acc_cor += ang_cos * (z_pos_eul * ang_cos - y_pos_eul * ang_sin)
        
        x_vel_cor = corr_velocity(delta_t, x_vel_eul, x_acc_eul)
        y_vel_cor = corr_velocity(delta_t, y_vel_eul, y_acc_eul)
        z_vel_cor = corr_velocity(delta_t, z_vel_eul, z_acc_eul)
        
        x_pos_cor = corr_position(delta_t, x_pos_eul, x_vel_eul)
        y_pos_cor = corr_position(delta_t, y_pos_eul, y_vel_eul)
        z_pos_cor = corr_position(delta_t, z_pos_eul, z_vel_eul)
        
        # Actual Acceleration, Velocity, and Position
        x_acc.append(pred_corr(x_acc_eul, x_acc_cor))
        y_acc.append(pred_corr(y_acc_eul, y_acc_cor))
        z_acc.append(pred_corr(z_acc_eul, z_acc_cor))
        
        x_vel.append(pred_corr(x_vel_eul, x_vel_cor))
        y_vel.append(pred_corr(y_vel_eul, y_vel_cor))
        z_vel.append(pred_corr(z_vel_eul, z_vel_cor))
        
        x_pos.append(pred_corr(x_pos_eul, x_pos_cor))
        y_pos.append(pred_corr(y_pos_eul, y_pos_cor))
        z_pos.append(pred_corr(z_pos_eul, z_pos_cor))

    return (time, x_acc, y_acc, z_acc, x_vel, y_vel, z_vel, 
            x_pos, y_pos, z_pos)    
#delta_t = z_pos[-1] * (time[-1] - time[-2]) / (z_pos[-2] - z_pos[-1])
#print("Interpolating: " + str(delta_t))
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Answer Questions
    
# Part A:     On Paper
    
# Part B:     Assume a projectile is launched due North at an initial height of
#             2 meters, level to the ground, and with an initial speed of 300
#             meters / second. What is the projectile's position when it hits
#             the ground? Neglect drag.
time,x_acc,y_acc,z_acc,x_vel,y_vel,z_vel,x_pos,y_pos,z_pos = projectile(
        height=2, initial_vel=300, north='North', drag=0)
print("North Projectile's Position:")
print("X Position: " + str(x_pos[-1]) + " m.")
print("Y Position: " + str(y_pos[-1]) + " m.")
print("Z Position: " + str(z_pos[-1]) + " m.")

plt.plot(x_pos, y_pos, 'g--')
plt.title("X - Y Path")
plt.xlabel("X Position (m)")
plt.ylabel("Y Position (m)")
# plt.savefig('Vertical Acceleration.pdf')
plt.show(block = True)

time,x_acc,y_acc,z_acc,x_vel,y_vel,z_vel,x_pos,y_pos,z_pos = projectile(
        height=2, initial_vel=300, north='East', drag=0)
print("\nEast Projectile's Position:")
print("X Position: " + str(x_pos[-1]) + " m.")
print("Y Position: " + str(y_pos[-1]) + " m.")
print("Z Position: " + str(z_pos[-1]) + " m.")

plt.plot(x_pos, y_pos, 'g--')
plt.title("X - Y Path")
plt.xlabel("X Position (m)")
plt.ylabel("Y Position (m)")
# plt.savefig('Vertical Acceleration.pdf')
plt.show(block = True)

    
time,x_acc,y_acc,z_acc,x_vel,y_vel,z_vel,x_pos,y_pos,z_pos = projectile()
    
# Part A:     Determine the x, y landing point of the pumpkin
#             Determined by when z was within 0.001 of 0
print("\nPumpkin Projectile:")
print("X Position: " + str(x_pos[-1])[:9] + " m.")
print("Y Position: " + str(y_pos[-1])[:8] + " m.")
print("Z Position: " + str(z_pos[-1])[:10] + " m.")

# Part B:     State altitude and Azimuth 
print("\nAltitude: " + str(altitude) + "° or " + str(alt))
print("Azimuth: " + str(azimuth) + "° or " + str(azi))

# Part C:     Flight time of the pumpkin to millisecond precision
print("\nTime: " + str(time[-1])[:6] + " s.")

# Part D:     Impact Velocity Vector to 0.01 m/s precision
print("\nX Velocity: " + str(x_vel[-1])[:7] + " m/s.")
print("Y Velocity: " + str(y_vel[-1])[:6] + " m/s.")
print("Z Velocity: " + str(z_vel[-1])[:8] + " m/s.")

total_vel = np.sqrt(x_vel[-1]**2 + y_vel[-1]**2 + z_vel[-1]**2)
print("Impact Velocity: " + str(total_vel)[:7] + " m/s.")

# Part E:     Create two plots. One to show the range vs height. The second
#             should show the vertical acceleration as a function of time.
#             Describe what is physically happening that produces the curve.

# Range Array 
r_arr = np.sqrt(np.array(x_pos)**2 + np.array(y_pos)**2)

# Plots the range vs height
plt.plot(r_arr, z_pos, 'g--')
plt.title("Range vs Height")
plt.xlabel("Range (m)")
plt.ylabel("Height (m)")
# plt.savefig('Range.pdf')
plt.show(block = True)

# Plots the vertical acceleration as a function of time
x_vel = [initial_vel * np.cos(alt) * np.sin(azi)]
y_vel = [initial_vel * np.cos(alt) * np.cos(azi)]

a = np.array(range(-130))

plt.plot(x_pos, y_pos, 'g--')
plt.plot(a*np.cos(azi)/np.sin(azi), a, 'b-')
plt.title("X - Y Path")
plt.xlabel("X Position (m)")
plt.ylabel("Y Position (m)")
# plt.savefig('Vertical Acceleration.pdf')
plt.show(block = True)
