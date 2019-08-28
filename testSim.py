# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:05:43 2019

@author: Aleks
"""

import numpy as np
import matplotlib.pyplot as plt
from plant import PneumCylinder, Valve


#%% System parameters
# VALVE
valve = Valve()
valve.loadRicherParam()

# CYLINDER
cyl = PneumCylinder()
cyl.loadRicherParam()


valve.Ps = 2e5

#%% Simulation parameters
Ts = 1e-4 # [s] step size
T_end = 0.2 # [s]
t = np.arange(0,T_end,Ts)

# Data containers
x = np.empty((4, len(t)+1))
xd = np.copy(x)

mdot = np.empty((4, len(t)))

### Initial condition
x[:,0] = np.array([0,0,cyl.Pa,cyl.Pa])

### Input
#u = 0*t # No input
#u = 5e-3 *np.sin(t) # Sinusoid
# Step input
u = np.ones(len(t))*(valve.Rh + valve.p)
t_step = 0.1
u[t < t_step] = u[t<t_step]*0
#u = -u

# Set to 0 or 1
plotValveArea1 = 0
plotValveArea2 = 0
fixPiston = 1


#%% Area test

"""
def AA1(x, R, p, nh):
    if x <= (p - R):
        A = 0
    elif x > (p - R) and x < (p + R):
        A = R**2 *np.arccos((p - x)/R) -(p - x)*np.sqrt(R**2 -(p - x)**2)
    elif x >= (p + R):
        A = np.pi*R**2
    else:
        raise ValueError('Something went wrong')
    return nh*A
def AA2(x, R, p, nh):
    xm = -x
    if xm <= (p - R):
        A = 0
    elif xm > (p - R) and xm < (p + R):
        A = R**2 *np.arccos((p - xm)/R) -(p - xm)*np.sqrt(R**2 -(p - xm)**2)
    elif xm >= (p + R):
        A = np.pi*R**2
    else:
        raise ValueError('Something went wrong')
    return nh*A
"""

# Show valve areas over the range
if plotValveArea1:
    xs = np.linspace(-1.1*(valve.p +valve.R), 1.1*(valve.p +valve.R), 40)
    Ao = np.empty(len(xs))
    Ao2 = np.empty(len(xs))
    for i in range(len(xs)):
        Ao[i] = valve.Avin(xs[i])
        Ao2[i] = valve.Avex(xs[i])
    
    plt.figure()
    plt.plot(xs*1e3,Ao*1e6,'bo-')
    plt.plot(xs*1e3,Ao2*1e6,'ro-')
    plt.ylabel('Area [mm$^2$]')
    plt.xlabel('Spool position [mm]')


# Show valve areas for the input sequence
if plotValveArea2:
    Ao = np.empty(len(u))
    Ao2 = np.empty(len(u))
    for i in range(len(u)):
        Ao[i] = valve.Avin(u[i])
        Ao2[i] = valve.Avex(u[i])
    
    plt.figure()
    plt.plot(t,Ao*1e6,'bo-')
    plt.plot(t,Ao2*1e6,'ro-')
    plt.ylabel('Area [mm^2]')
    plt.xlabel('Spool position [mm]')


#%% Simulation


for i in range(len(t)):
#for i in range(10):
    # Update mass flow rates
    
    mdot[0,i] = valve.mdot(valve.Ps, x[2,i], valve.Avin(u[i])) # Supply -> p1
    mdot[1,i] = valve.mdot(x[2,i], valve.Pa, valve.Avex(u[i])) # p1 -> Exhaust
    mdot[2,i] = valve.mdot(valve.Ps, x[3,i], valve.Avex(u[i])) # Supply -> p2
    mdot[3,i] = valve.mdot(x[3,i], valve.Pa, valve.Avin(u[i])) # p2 -> Exhaust
    
    # Calculate derivative
    xd[:,i] = cyl.derivative(x[:,i], mdot[0,i], mdot[1,i], mdot[2,i], mdot[3,i])
    
    if fixPiston: # Fixate piston
        xd[0,i] = 0
        xd[1,i] = 0
    else:
        ### Endstop test
        print(x[0,i])
        if np.abs(x[0,i]) >= cyl.S and x[0,i]*xd[0,i] > 0:
            #
            # If the piston is at either end and going further away
            xd[0,i] = 0 # No change in position
            x[1,i] = 0 # Set velocity to zero (abrupt stop)
        
    # Integrate
    x[:,i+1] = x[:,i] + Ts*xd[:,i]

# Remove last entries of x to match the time array
x = x[:,0:-1]
    

#%%
fig,ax = plt.subplots(5,1, sharex=True, figsize=(12,10))
ax[0].plot(t, x[0]*1e3, label='x [mm]')
ax[0].plot(t, x[1]*1e3, ls='--', label='xd [mm/s]')
ax[0].set_ylabel('Motion')
ax[0].legend()

ax[1].plot(t, x[2]*1e-5, label='P1')
ax[1].plot(t, x[3]*1e-5, '--',label='P2')
ax[1].set_ylabel('[bar]')
ax[1].legend()

ax[2].plot(t, u, label='Spool pos [mm]')
ax[2].set_ylabel('Input [mm]')
ax[2].set_ylim([-1.1*(valve.Rh+valve.p), 1.1*(valve.Rh+valve.p)])
ax[2].legend()

ax[3].plot(t, mdot[0]*1e3, color='g', label='In')
ax[3].plot(t, mdot[1]*1e3, ls='--', color='orange', label='Out')
ax[3].legend()
ax[3].set_ylabel('$\dot{m}_1$ [g/s]')

ax[4].plot(t, mdot[2]*1e3, color='g')
ax[4].plot(t, mdot[3]*1e3, ls='--', color='orange')
ax[4].set_ylabel('$\dot{m}_2$ [g/s]')
ax[4].set_xlabel('Time [s]')