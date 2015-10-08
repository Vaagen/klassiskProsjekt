import numpy as np
import matplotlib
import matplotlib.pyplot as plt


from math import *    #this command gives you acces to math functions, such as sin(), pow() etc



G = 1    # Gravitational constant
M = 10    # Central Mass
m1 = 0.0001    # mass of particle 1
m2 = 0.0001    # mass of particle 2
dt = 0.00001  # integration timestep

# create files to save simulation data in:
f = open('tertiarypos1.txt','w') # notice: the write option 'w' erases previous data in the file
f2 = open('tertiarypos2.txt','w')
f3 = open('tertiaryenergy1.txt','w')
f4 = open('tertiaryenergy2.txt','w')
f5 = open('tertiaryangularmomentum1.txt','w')
f6 = open('tertiaryangularmomentum2.txt','w')

x1 = 1             #initial position in x-direction, particle 1
y1 = 0             #initial position in y-direction, particle 1
vx1 = 0           #initial velocity in x-direction, particle 1
vy1 = 1            #initial velocity in y-direction, particle 1

x2 = -1             #initial position in x-direction, particle 2
y2 = 0             #initial position in y-direction, particle 2
vx2 = 0            #initial velocity in x-direction, particle 2
vy2 = 1            #initial velocity in y-direction, particle 2

time=0.0          #this is the start time
endtime=200000*dt        #total simulation time

i = 1 # variable used to save only each <plotdensity>'th value...
plotdensity = 1000


def r1():    # distance between origin and particle 1
    return (x1**2 + y1**2)**0.5
def r2():    # distance between origin and particle 2
    return (x2**2 + y2**2)**0.5
def r1r2():  # distance between particle 1 and 2
    return ( (x1 - x2)**2 + (y1 - y2)**2)**0.5

def gradV1():   # part of gradV1 in common for both x- and y-direction
    return (G*M*m1/r1()**(3/2) + G*m1*m2/r1r2()**(3/2))
def gradV1x():  # x component of gradV1
    return gradV1()*x1
def gradV1y():  # y component of gradV1
    return gradV1()*y1

def gradV2():   # part of gradV2 in common for both x- and y-direction
    return (G*M*m2/r2()**(3/2) + G*m1*m2/r1r2()**(3/2))
def gradV2x(): # x component of gradV2
    return gradV1()*x2
def gradV2y():  # y component of gradV2
    return gradV1()*y2

while (time < endtime):

        time = time + dt
        i = i+1

        fx1 = -gradV1x()/m1
        fy1 = -gradV1y()/m1
        vxm1 = vx1 + 0.5*dt*fx1
        vym1 = vy1 + 0.5*dt*fy1

        fx2 = -gradV2x()/m2
        fy2 = -gradV2y()/m2
        vxm2 = vx2 + 0.5*dt*fx2
        vym2 = vy2 + 0.5*dt*fy2

        x1 = x1 + vxm1*dt
        y1 = y1 + vym1*dt
        fx1 = -gradV1x()/m1
        fy1 = -gradV1y()/m1
        vx1 = vxm1 + 0.5*dt*fx1
        vy1 = vym1 + 0.5*dt*fy1

        x2 = x2 + vxm2*dt
        y2 = y2 + vym2*dt
        fx2 = -gradV2x()/m2
        fy2 = -gradV2y()/m2
        vx2 = vxm2 + 0.5*dt*fx2
        vy2 = vym2 + 0.5*dt*fy2

        if (i % plotdensity == 0):
            E1 = 0.5*m1*(vx1**2 + vy1**2) - G*M*m1/r1() - G*m1*m2/r1r2() # total energy of particle 1
            E2 = 0.5*m2*(vx2**2 + vy2**2) - G*M*m2/r2() - G*m1*m2/r1r2() # total energy of particle 2
            
            L1 = m1*(x1*vy1 - y1*vx1)   # total angular momentum of particle 1, in z-direction
            L2 = m2*(x2*vy2 - y2*vx2)   # total angular momentum of particle 2, in z-direction
            
            f.write("%f %f\n" % (x1,y1))
            f2.write("%f %f\n" % (x2,y2))
            f3.write("%f %f\n" % (time,E1))
            f4.write("%f %f\n" % (time, E2))
            f5.write("%f %f\n" % (time, L1))
            f6.write("%f %f\n" % (time, L2))

#closing files
f.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()



        