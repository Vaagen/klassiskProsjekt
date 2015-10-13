import numpy as np
import matplotlib
import matplotlib.pyplot as plt


from math import *    #this command gives you acces to math functions, such as sin(), pow() etc

        
k=1    # force constant
m=1.0    # particle mass relative to particle mass

X=100000             #initial position in x-direction,
Y=0             #initial position in y-direction
Vx=0            #initial velocity in x-direction
Vy=20000            #initial velocity in y-direction

dT=0.01  # integration timestep
Time=0.0          #this is the start time
endTime=650*dT        #total simulation time

r0 = (X**2+Y**2)**0.5      #rewriting to dimensionless parameters
v0 = (Vx**2 + Vy**2)**0.5
x = X /r0
y = Y / r0
vx = Vx / v0
vy = Vy / v0
dt = v0/r0*dT
time = v0/r0*Time
endtime = v0/r0*endTime

# create files to save simulation data in:
f = open('binarypos.txt','w') # notice: the write option 'w' erases previous data in the file
f2 = open('binaryenergy.txt','w')
f3 = open('binaryangularmomentum.txt','w')

maxEnergyDeviation = 0              # maximum deviation of energy for a particle, relative to start
maxAngularMomentumDeviation = 0     # maximum deviation of angular momentum for a particle, relative to start

startEnergy=0.5*m*(vx**2+vy**2) - k/(x**2+y**2)**(0.5)
startAngMomentum = m*(x*vy - y*vx)

while (time < endtime):  

           time=time+dt

           fx=-k*x/(m*(x**2+y**2)**(3/2) )          # calculate force/mass on the particle at time "t" in x-direction
           fy=-k*y/(m*(x**2+y**2)**(3/2) )          # calculate force/mass on the particle at time "t" in y-direction
           vxm=vx+0.5*dt*fx    # first step in Verlet algorithm : calculate velocity at time  "t+dt/2", x-direction
           vym=vy+0.5*dt*fy    # first step in Verlet algorithm : calculate velocity at time  "t+dt/2", y-direction

           x=x+vxm*dt          # calculate new position x at time "t+dt", x-direction
           y=y+vym*dt          # calculate new position x at time "t+dt", y-direction
           fx=-k*x/(m*(x**2+y**2)**(3/2) )          # calculate force using the new postion x at "t+dt", x-direction
           fy=-k*y/(m*(x**2+y**2)**(3/2) )          # calculate force using the new postion x at "t+dt", y-direction
           vx=vxm+0.5*dt*fx    # calculate velocity at "t+dt", x-direction
           vy=vym+0.5*dt*fy    # calculate velocity at "t+dt", y-direction


           energy=0.5*m*(vx**2+vy**2) - k/(x**2+y**2)**(0.5)    # calculate the total energy, to check that it is approximately conserved!
           angMomentum = m*(x*vy - y*vx)        #calculate angular momentum in z-direction
          

 
           f.write("%f %f\n" % (x, y))            # Writes x and y positions to file 1       
           f2.write("%f %f\n" % (time, energy))   #  write total energy  to file 2 (can be used to plot total energy as a  function of time)
           f3.write("%f %f\n" % (time, angMomentum)) # writes total angular momentum to file 3

           if (abs(1-energy/startEnergy) > maxEnergyDeviation):
                    maxEnergyDeviation = 1 - energy/startEnergy
           if (abs(1-angMomentum/startAngMomentum) > maxAngularMomentumDeviation):
                    maxAngularMomentumDeviation = angMomentum/startAngMomentum - 1

print("Maximum deviation of energy for a particle relative to start energy is %f\n" % maxEnergyDeviation)
print("Maximum deviation of angular momentum for a particle relative to start angular momentum is %f\n" % maxAngularMomentumDeviation)

# close files
f.close()
f2.close()
f3.close()
