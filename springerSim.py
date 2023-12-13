#!/usr/bin/env python3

import numpy as np
from math import ceil
from time import perf_counter
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Settings

barrelDiameter_m = 0.0133 # m

dartMass_kg = 0.0012 # kg
dartStaticFriction_N = 2 # N
dartDynamicFriction_N = 1 # N

dt_s = 0.00001 # s
maxTime_s = 0.05 # s

temperature_C = 25 # C

# Caliburn
barrelLength_m = 0.33 # 13in
plungerDiameter_m = 0.034798 # m
plungerDraw_m = 0.15138
precompression_m = 0.040818
#springRate_Npm = 679.651333723 # N/m k26
springForce_N = 13.3*9.8 # k26 N at starting distance
#springForce_N = 6.51*9.8 # 788 N at starting distance
springMass_kg = 0.0284
plungerRodMass_kg = 0.03 # kg
plungerFriction_N = 1 # N
plungerDeadspace_m3 = 0.000004 # cubic meters
barrelDeadspace_m3 = 0.000004 # cubic meters


'''
# Pistol
barrelLength_m = 0.1016 # m 
plungerDiameter_m = 0.03495 # m
plungerDraw_m = 0.0508
precompression_m = 0
#springRate_Npm = 679.651333723 # N/m k26
#springForce_N = 24*9.8 # K14
springMass_kg = 0.0284
plungerRodMass_kg = 0.04 # kg
plungerFriction_N = 0.5 # N
'''

# Sentinel
#plungerRodMass_kg = 0.03
#plungerDiameter_m = 1.25*.0254
#plungerDraw_m = 3.25*.0254

# Nite Finder
#plungerRodMass_kg = 0.0001
#plungerDiameter_m = 1*.0254
#plungerDraw_m = 2*.0254
#precompression_m = 0
#deadspace = 0
#barrelLength_m = 12*.0254
#springRate_Npm = 2*9.8/plungerDraw_m

pi = 3.14159265358979323
pa = 101325
ft = 3.28084


def airDensity_kgpm3(temperature_C, pressure):
    R=287.05 # Specific gas constant, J/(kg*K)
    C_K_offset=273.15
#    print((pressure+1)*pa, "Pa", temperature_C+C_K_offset,"K")
    return (pressure+1)*pa/(R*(temperature_C+C_K_offset))

def simulate():
    arrayLength = ceil(maxTime_s/dt_s)+2
    times_s = np.full(arrayLength, np.nan, dtype=float)
    averagePressures_atm = np.full(arrayLength, np.nan, dtype=float)
    plungerPressures_atm = np.full(arrayLength, np.nan, dtype=float)
    barrelPressures_atm = np.full(arrayLength, np.nan, dtype=float)
    plungerVolumes_m3 = np.full(arrayLength, np.nan, dtype=float)
    barrelVolumes_m3 = np.full(arrayLength, np.nan, dtype=float)
    plungerPositions_m = np.full(arrayLength, np.nan, dtype=float)
    plungerVelocities_mps = np.full(arrayLength, np.nan, dtype=float)
    dartPositions_m = np.full(arrayLength, np.nan, dtype=float)
    dartVelocities_mps = np.full(arrayLength, np.nan, dtype=float)
    bernoulliPressureDifferences_atm = np.full(arrayLength, np.nan, dtype=float)

    PE_Spring_J = np.full(arrayLength, np.nan, dtype=float)
    KE_Plunger_J = np.full(arrayLength, np.nan, dtype=float)
    KE_Dart_J = np.full(arrayLength, np.nan, dtype=float)
    KE_PlungerTubeAir_J = np.full(arrayLength, np.nan, dtype=float)
    KE_BarrelAir_J = np.full(arrayLength, np.nan, dtype=float)
    PE_PlungerTubePressure_J = np.full(arrayLength, np.nan, dtype=float)
    PE_BarrelPressure_J = np.full(arrayLength, np.nan, dtype=float)
    TE_PlungerFriction_J = np.full(arrayLength, np.nan, dtype=float)
    TE_DartFriction_J = np.full(arrayLength, np.nan, dtype=float)
    TE_AirFriction_J = np.full(arrayLength, np.nan, dtype=float)
    totalEnergy_J = np.full(arrayLength, np.nan, dtype=float)

    starttime = perf_counter()
    plungerRadius_m = plungerDiameter_m/2
    barrelRadius_m = barrelDiameter_m/2
    plungerStart_m = precompression_m + plungerDraw_m
    effectivePlungerMass_kg = plungerRodMass_kg + springMass_kg/2
    springRate_Npm = springForce_N/plungerStart_m # N/m
    print(plungerStart_m, springRate_Npm)
#    plungerForce_N = plungerStart_m * springRate_Npm - plungerFriction_N # compression spring
    plungerForce_N = springForce_N - plungerFriction_N # constant force spring
    springPotentialEnergy_J = 0.5*springRate_Npm*plungerStart_m**2-0.5*springRate_Npm*precompression_m**2
    
    times_s[0] = 0
    plungerPressures_atm[0] = 0
    barrelPressures_atm[0] = 0
    averagePressures_atm[0] = 0
    plungerVolumes_m3[0] = plungerDeadspace_m3 + (plungerStart_m-precompression_m) * pi * plungerRadius_m**2
    barrelVolumes_m3[0] = barrelDeadspace_m3
    dartPositions_m[0] = 0
    dartVelocities_mps[0] = 0
    plungerPositions_m[0] = plungerStart_m
    plungerVelocities_mps[0] = 0
    bernoulliPressureDifferences_atm[0] = 0

    i=0
    airVolume = plungerVolumes_m3[0]+barrelVolumes_m3[0]
    
    while dartPositions_m[i] < barrelLength_m and times_s[i] < maxTime_s and (dartVelocities_mps[i] > 0 or dartPositions_m[i] == 0):
        i += 1
        times_s[i] = times_s[i-1]+dt_s
    
        
        if plungerPositions_m[i-1] > precompression_m:
            plungerVelocities_mps[i] = plungerVelocities_mps[i-1] + plungerForce_N/effectivePlungerMass_kg*dt_s
            plungerPositions_m[i] = max(precompression_m, plungerPositions_m[i-1] - plungerVelocities_mps[i]*dt_s)
            if plungerPositions_m[i] == precompression_m:
                print("plunger hit end at", dartPositions_m[i-1]*ft*12)
        else:
            plungerVelocities_mps[i] = 0
            plungerPositions_m[i] = precompression_m
        

        plungerVolumes_m3[i] = plungerDeadspace_m3 + (plungerPositions_m[i]-precompression_m) * pi * plungerRadius_m**2
        barrelVolumes_m3[i] = barrelDeadspace_m3 + dartPositions_m[i-1] * pi * barrelRadius_m**2
        averagePressures_atm[i] = airVolume/(plungerVolumes_m3[i]+barrelVolumes_m3[i])-1

        ptAirDensity_kgpm3 = airDensity_kgpm3(temperature_C, plungerPressures_atm[i-1])
        barrelAirDensity_kgpm3 = airDensity_kgpm3(temperature_C, barrelPressures_atm[i-1])
        bernoulliPressureDifferences_atm[i] = (0.5 * ptAirDensity_kgpm3 * plungerVelocities_mps[i-1]**2 - 0.5 * barrelAirDensity_kgpm3 * (plungerVelocities_mps[i-1]*plungerRadius_m**2/barrelRadius_m**2)**2)/pa
#        bernoulliPressureDifferences_atm[i] = (0.5 * airDensity_kgpm3(temperature_C, averagePressures_atm[i-1]) * plungerVelocities_mps[i-1]**2 - 0.5 * airDensity_kgpm3(temperature_C, averagePressures_atm[i-1]) * (plungerVelocities_mps[i-1]*plungerRadius_m**2/barrelRadius_m**2)**2)/pa
#        bernoulliPressureDifferences_atm[i] = 0
        # averagePressure = (plungerVolume * plungerPressure + barrelVolume * barrelPressure) / (plungerVolume+barrelVolume)

        plungerPressures_atm[i] = averagePressures_atm[i] + barrelVolumes_m3[i]*bernoulliPressureDifferences_atm[i]/(plungerVolumes_m3[i]+barrelVolumes_m3[i])
        
        barrelPressures_atm[i] = averagePressures_atm[i] - plungerVolumes_m3[i]*bernoulliPressureDifferences_atm[i]/(plungerVolumes_m3[i]+barrelVolumes_m3[i])

#        print("d", round(bernoulliPressureDifferences_atm[i],2), "a", round(averagePressures_atm[i],2), "x", round(plungerPressures_atm[i],2), "y", round(barrelPressures_atm[i],2), plungerVelocities_mps[i], dartVelocities_mps[i-1])
#        print("a", round(averagePressures_atm[i],2), "b", plungerVolumes_m3[i], "c", barrelVolumes_m3[i], "d", round(bernoulliPressureDifferences_atm[i],3), "x", round(plungerPressures_atm[i],2), "y", round(barrelPressures_atm[i],2))

        if plungerVelocities_mps[i] >= 0:
            plungerForce_N = springForce_N \
                - plungerFriction_N \
                - plungerPressures_atm[i] * pa * pi * plungerRadius_m**2
        else:
            plungerForce_N = springForce_N \
                + plungerFriction_N \
                - plungerPressures_atm[i] * pa * pi * plungerRadius_m**2
            
        if dartVelocities_mps[i-1] > 0:
            dartForce = barrelPressures_atm[i] * pa * pi * barrelRadius_m**2 - dartDynamicFriction_N
        else:
            dartForce = max(0, barrelPressures_atm[i] * pa * pi * barrelRadius_m**2 - dartStaticFriction_N)
            
        dartVelocities_mps[i] = dartVelocities_mps[i-1] + dartForce/dartMass_kg*dt_s
        dartPositions_m[i] = dartPositions_m[i-1] + dartVelocities_mps[i]*dt_s

        #Potential Energy
        PE_Spring_J[i] = 0.5*springRate_Npm*plungerPositions_m[i]**2-0.5*springRate_Npm*precompression_m**2
        PE_PlungerTubePressure_J[i] = plungerPressures_atm[i]*pa*plungerVolumes_m3[i]
        PE_BarrelPressure_J[i] = barrelPressures_atm[i]*pa*barrelVolumes_m3[i]

        #Kinetic Energy
        KE_Plunger_J[i] = 0.5*effectivePlungerMass_kg*plungerVelocities_mps[i]**2
        KE_Dart_J[i] = 0.5*dartMass_kg*dartVelocities_mps[i]**2
        KE_PlungerTubeAir_J[i] = 0.5 * ptAirDensity_kgpm3*plungerVolumes_m3[i] * plungerVelocities_mps[i]**2
        KE_BarrelAir_J[i] = 0.5 * barrelAirDensity_kgpm3*barrelVolumes_m3[i] * dartVelocities_mps[i]**2

        #Thermal Energy
        TE_PlungerFriction_J[i] = plungerFriction_N*(plungerStart_m-plungerPositions_m[i])
        TE_DartFriction_J[i] = dartDynamicFriction_N*dartPositions_m[i]
        TE_AirFriction_J[i] = 0

        totalEnergy_J[i] = PE_Spring_J[i] + PE_PlungerTubePressure_J[i] + PE_BarrelPressure_J[i] + KE_Plunger_J[i] + KE_Dart_J[i] + KE_PlungerTubeAir_J[i] + KE_BarrelAir_J[i] + TE_PlungerFriction_J[i] + TE_DartFriction_J[i] + TE_AirFriction_J[i]

        
#        print(plungerPositions_m[i]-plungerStart_m, plungerVolumes_m3[i], plungerPressures_atm[i], barrelPressures_atm[i], dartVelocities_mps[i])
        
    print(i, "cycles, max", arrayLength)
    print(perf_counter()-starttime, "s runtime")
    print(round(dartVelocities_mps[i]*ft,3), "fps in", round(times_s[i],2), "s and", round(dartPositions_m[i]*ft*12,1),"in")
    print(np.amax(plungerPressures_atm),np.amax(barrelPressures_atm),np.amax(averagePressures_atm),np.amax(bernoulliPressureDifferences_atm))
    
    maxFPS[j] = np.amax(dartVelocities_mps[0:i])*ft
#    maxP[j] = np.amax(pressures[0:i])*pa/1000
    idealBarrelLength[j] = dartPositions_m[np.argmax(dartVelocities_mps[0:i])]
    try:
        plungerStop_t = np.nonzero(plungerPositions_m <= precompression_m)[0][0]
    except IndexError:
        plungerStop_t = len(plungerPositions_m)
    print(round(maxFPS[j],3), "fps peak at", round(idealBarrelLength[j]*1000), "mm")
    print(round(springPotentialEnergy_J,2), "J", round(totalEnergy_J[i],2), "J")
    print("plunger stops at", plungerStop_t)

    plt.ioff()
    fig = plt.figure(figsize=(18,12), dpi=600)
    plt.title('SpringerSim')
    plt.ylabel('fps / mm / kPa')
#    variableName = "Plunger Draw: "+str(round(v*1000,1))+" mm "
    variableName = ""
#    xAxis = times_s
#    plt.xlabel('Time elapsed (seconds)')
    xAxis = dartPositions_m*ft*12
    plt.xlabel('Dart travel (inches)')
    plt.plot(xAxis, dartVelocities_mps*ft, label=variableName+"velocity, max "+str(round(maxFPS[j],1))+"fps @ "+str(round(idealBarrelLength[j]*1000))+" mm")
    plt.plot(xAxis[:plungerStop_t], (plungerPositions_m[:plungerStop_t]-precompression_m)*1000, label=variableName+"plunger position (mm)")
#    plt.plot(xAxis, dartPositions_m*200, label=variableName+"dart position (mm)")
    plt.plot(xAxis, barrelPressures_atm*101.325, label=variableName+"barrel pressure (kPa)")
    plt.plot(xAxis, plungerPressures_atm*101.325, label=variableName+"plunger pressure (kPa)")
    
    plt.legend(loc='center right', markerscale=100)
    fig.savefig('springerSim '+str(datetime.now().replace(microsecond=0)).replace(':','-')+'.png', bbox_inches='tight')

    plt.clf()

    plt.ylabel('Joules')
    plt.xlabel('Time elapsed (seconds)')
    plt.stackplot(xAxis, PE_Spring_J, PE_PlungerTubePressure_J, PE_BarrelPressure_J, KE_PlungerTubeAir_J, KE_BarrelAir_J, TE_PlungerFriction_J, TE_DartFriction_J, TE_AirFriction_J, KE_Plunger_J, KE_Dart_J,
        labels=['PE_Spring', 'PE_PlungerTubePressure', 'PE_BarrelPressure', 'KE_PlungerTubeAir', 'KE_BarrelAir', 'TE_PlungerFriction', 'TE_DartFriction', 'TE_AirFriction', 'KE_Plunger', 'KE_Dart'])
    plt.legend(loc='lower right', markerscale=100)
    plt.ylim(ymin=0)
    fig.savefig('springerEnergyStacked '+str(datetime.now().replace(microsecond=0)).replace(':','-')+'.png', bbox_inches='tight')

    plt.clf()
    plt.ylabel('Joules')
    plt.xlabel('Time elapsed (seconds)')
    plt.plot(xAxis, PE_Spring_J, label="PE_Spring")
    plt.plot(xAxis, PE_PlungerTubePressure_J, label="PE_PlungerTubePressure")
#    plt.plot(xAxis, PE_BarrelPressure_J, label="PE_BarrelPressure") 
#    plt.plot(xAxis, PE_PlungerTubePressure_J+PE_BarrelPressure_J, label="PE_TotalPressure")
#    plt.plot(xAxis, KE_PlungerTubeAir_J, label="KE_PlungerTubeAir") 
#    plt.plot(xAxis, KE_BarrelAir_J) 
#    plt.plot(xAxis, TE_PlungerFriction_J, label="TE_PlungerFriction")
#    plt.plot(xAxis, TE_DartFriction_J) 
#    plt.plot(xAxis, TE_AirFriction_J) 
    plt.plot(xAxis, KE_Plunger_J, label="KE_Plunger") 
#    plt.plot(xAxis, KE_Dart_J)
    plt.plot(xAxis, totalEnergy_J-springPotentialEnergy_J, label="Error")
    plt.legend(loc='lower right', markerscale=100)
    plt.ylim(ymin=0)
    fig.savefig('springerEnergyLines '+str(datetime.now().replace(microsecond=0)).replace(':','-')+'.png', bbox_inches='tight')


iterationArray1 = np.arange(plungerDraw_m/2,plungerDraw_m,0.01)
maxFPS = np.zeros(len(iterationArray1), dtype=float)
maxP = np.zeros(len(iterationArray1), dtype=float)
idealBarrelLength = np.zeros(len(iterationArray1), dtype=float)
for j, v in enumerate(iterationArray1):
    plungerDraw_m = v
#    simulate()
simulate()




'''
plt.plot(iterationArray1*1000, maxFPS, label="max FPS")
plt.plot(iterationArray1*1000, idealBarrelLength*1000, label="ideal barrel length, mm")
#plt.plot(iterationArray1*1000, maxP*10, label="max pressure")
plt.title('SpringerSim2')
plt.ylabel('fps / mm / kPa')
plt.xlabel('plunger draw (mm)')
plt.legend(loc='lower right', markerscale=100)
#plt.ylim(ymin=0)
fig.savefig('springerSimMetas.png', bbox_inches='tight')
'''
#print(idealBarrelLength/iterationArray1)