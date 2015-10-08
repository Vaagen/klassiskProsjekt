import numpy as np
import matplotlib
import matplotlib.pyplot as plt


from math import *    #this command gives you acces to math functions, such as sin(), pow() etc



k=1    # spring constant, equal to G*M
m1 = 1.0    # mass of particle 1
m2 = 1.0    # mass of particle 2
dt=0.02  # integration timestep

# create files to save simulation data in:
f = open('tertiarypos1.txt','w') # notice: the write option 'w' erases previous data in the file
f2 = open('tertiarypos2.txt','w')
f3 = open('tertiaryenergy1.txt','w')
f4 = open('tertiaryenergy2.txt','w')

x1 = 1             #initial position in x-direction, particle 1
y1 = 0             #initial position in y-direction, particle 1
vx1 = 0            #initial velocity in x-direction, particle 1
vy1 = 1            #initial velocity in y-direction, particle 1

x2 = 0             #initial position in x-direction, particle 2
y2 = 1             #initial position in y-direction, particle 2
vx2 = -1            #initial velocity in x-direction, particle 2
vy2 = 0            #initial velocity in y-direction, particle 2

time=0.0          #this is the start time
endtime=350*dt        #total simulation time


while (time < endtime)

        time = time + dt

        