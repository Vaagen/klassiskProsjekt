import numpy as np
import matplotlib
import matplotlib.pyplot as plt


from math import *    #this command gives you acces to math functions, such as sin(), pow() etc



k=1    # spring constant, equal to G*M
m1 = 1.0    # mass of particle 1
m2 = 1.0    # mass of particle 2
dt=0.02  # integration timestep

# create files to save simulation data in:
f = open('xyposition.txt','w') # notice: the write option 'w' erases previous data in the file
f2 = open('energy.txt','w')