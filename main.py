import numpy as np
import matplotlib
import matplotlib.pyplot as plt


from math import *    #this command gives you acces to math functions, such as sin(), pow() etc

        
k=1    # spring constant
m=1.0    # particle mass 
dt=0.002  # integration timestep

# create files to save simulation data in:
f = open('xyposition.txt','w') # notice: the write option 'w' erases previous data in the file 
f2 = open('energy.txt','w')



x=1             #initial position in x-direction
y=0             #initial position in y-direction
vx=0            #initial velocity in x-direction
vy=1            #initial velocity in y-direction
time=0.0          #this is the start time
endtime=3000*dt        #total simulation time



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
          

 
           f.write("%f %f\n" % (x, y))            # Writes x and y positions to file 1       
           f2.write("%f %f\n" % (time, energy))   #  write total energy  to file 2 (can be used to plot total energy as a  function of time)

# close files 
f.close()
f2.close()
