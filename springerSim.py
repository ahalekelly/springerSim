#!/usr/bin/env python3

import numpy as np
from math import ceil
from time import clock
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.ioff()
fig = plt.figure(figsize=(18,12), dpi=600)
plt.title('Springer Simulation')
plt.ylabel('fps / mm / MPa')

# Settings

# Caliburn
barrelLength = 0.321754 # m 
barrelLength = 3.5
plungerDiameter = 0.03495 # m
barrelDiameter = 0.0127 # m
plungerDraw = 0.15138
precompression = 0.040818
#springRate = 679.651333723 # N/m k26
springForce = 13.3*9.8 # N at starting distance
springMass = 0.0284
plungerRodMass = 0.04 # kg
plungerFriction = 0.5 # N
deadspace = 0.00000792152 # cubic meters

# Sentinel
#plungerMass = 0.03
plungerDiameter = 1.25*.0254
plungerDraw = 3.25*.0254

# Nite Finder
#plungerMass = 0.0001
#plungerDiameter = 1*.0254
#plungerDraw = 2*.0254
#precompression = 0
#deadspace = 0
#barrelLength = 12*.0254
#springRate = 2*9.8/plungerDraw

dartMass = 0.0012 # kg
dartFriction = 2.5 # N

dt = 0.00001
maxTime = 0.1


pi = 3.14159265358979323
atm = 101325
ft = 3.28084

arrayLength = ceil(maxTime/dt)+2
times = np.full(arrayLength, np.nan, dtype=float)
pressures = np.full(arrayLength, np.nan, dtype=float)
volumes = np.full(arrayLength, np.nan, dtype=float)
plungerPositions = np.full(arrayLength, np.nan, dtype=float)
plungerVelocities = np.full(arrayLength, np.nan, dtype=float)
dartPositions = np.full(arrayLength, np.nan, dtype=float)
dartVelocities = np.full(arrayLength, np.nan, dtype=float)

def simulate():
    timestart = clock()
    plungerRadius = plungerDiameter/2
    barrelRadius = barrelDiameter/2
    plungerStart = precompression + plungerDraw
    plungerMass = plungerRodMass + springMass/2
    springRate = springForce/plungerStart # N/m
    print(plungerStart, springRate)
    plungerForce = plungerStart * springRate - plungerFriction
    airVolume = deadspace + (plungerStart-precompression) * pi * plungerRadius**2    
    
    times[0] = 0
    pressures[0] = 0
    volumes[0] = airVolume
    dartPositions[0] = 0
    dartVelocities[0] = 0
    plungerPositions[0] = plungerStart
    plungerVelocities[0] = 0
    i=0
    
    while dartPositions[i] < barrelLength and times[i] < maxTime and (dartVelocities[i] > 0 or dartPositions[i] == 0):
        i += 1
        times[i] = times[i-1]+dt
    
        
        if plungerPositions[i-1] > precompression:
            plungerVelocities[i] = plungerVelocities[i-1] + plungerForce/plungerMass*dt
            plungerPositions[i] = max(precompression, plungerPositions[i-1] - plungerVelocities[i]*dt)
            if plungerPositions[i] == precompression:
                print("plunger hit end at", dartPositions[i-1]*ft*12)
        else:
            plungerVelocities[i] = 0
            plungerPositions[i] = precompression
        
        volumes[i] = deadspace \
            + (plungerPositions[i]-precompression) * pi * plungerRadius**2 \
            + dartPositions[i-1] * pi * barrelRadius**2
        
        pressures[i] = airVolume/volumes[i]-1
        if plungerVelocities[0] >= 0:
            plungerForce = plungerPositions[i] * springRate \
                - plungerFriction \
                - pressures[i] * atm * pi * plungerRadius**2
        else:
            plungerForce = plungerPositions[i] * springRate \
                + plungerFriction \
                - pressures[i] * atm * pi * plungerRadius**2
            
        
        if dartVelocities[i-1] > 0:
            dartForce = pressures[i] * atm * pi * barrelRadius**2 - dartFriction
        else:
            dartForce = max(0, pressures[i] * atm * pi * barrelRadius**2 - dartFriction)
            
        dartVelocities[i] = dartVelocities[i-1] + dartForce/dartMass*dt
        dartPositions[i] = dartPositions[i-1] + dartVelocities[i]*dt
        
        
#        print(i, plungerPositions[i], dartForce, dartVelocities[i], dartPositions[i])
        
    print(i, "/", arrayLength)
    print(clock()-timestart, "s runtime")
    print(dartVelocities[i]*ft, "fps in", times[i], "s")

    xAxis = times
    plt.xlabel('Time elapsed (seconds)')
    xAxis = dartPositions*ft*12
    plt.xlabel('Dart travel (inches)')
    plt.plot(xAxis, dartVelocities*ft, label="velocity, max "+str(int(np.amax(dartVelocities[0:i])*ft))+"fps @ "+str(round(dartPositions[np.argmax(dartVelocities[0:i])]*ft*12, 1))+" in")
    plt.plot(xAxis, (plungerPositions-precompression)*1000, label="plunger position (mm)")
#    plt.plot(xAxis, dartPositions*200, label="plunger position (mm)")
    plt.plot(xAxis, pressures*atm/1000, label="pressure (MPa)")
    
springForce=8*9.8
simulate()
springForce=10*9.8
simulate()
springForce=13*9.8
simulate()

#plt.ylim(0,80)

plt.legend(loc='lower right', markerscale=10)
fig.savefig('springerSim.png', bbox_inches='tight')