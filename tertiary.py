import numpy as np
import matplotlib
import matplotlib.pyplot as plt


from math import *    #this command gives you acces to math functions, such as sin(), pow() etc



G = 1          # Gravitational constant
M = 2         # Central Mass
m1 = 0.000001    # mass of particle 1
m2 = 0.0    # mass of particle 2
dt = 0.0001   # integration timestep

x1 = 1             #initial position in x-direction, particle 1
y1 =0              #initial position in y-direction, particle 1
vx1 = 0            #initial velocity in x-direction, particle 1
vy1 = 0.5            #initial velocity in y-direction, particle 1

x2 = 1            #initial position in x-direction, particle 2
y2 = 1             #initial position in y-direction, particle 2
vx2 = 0            #initial velocity in x-direction, particle 2
vy2 = 0            #initial velocity in y-direction, particle 2

time=0.0                 #this is the start time
endtime=20        #total simulation time

i = 1                # variable used to save only each <plotdensity>'th value...
plotdensity = 10

maxEnergyDeviation = 0              # maximum deviation of energy for a particle, relative to start
maxAngularMomentumDeviation = 0     # maximum deviation of angular momentum for a particle, relative to start


# create files to save simulation data in:
f = open('tertiarypos1.txt','w') # notice: the write option 'w' erases previous data in the file
f2 = open('tertiarypos2.txt','w')
f3 = open('tertiaryenergy1.txt','w')
f4 = open('tertiaryenergy2.txt','w')
f5 = open('tertiaryangularmomentum1.txt','w')
f6 = open('tertiaryangularmomentum2.txt','w')



def r1():    # distance between origin and particle 1
    return (x1**2 + y1**2)**0.5
def r2():    # distance between origin and particle 2
    return (x2**2 + y2**2)**0.5
def r1r2():  # distance between particle 1 and 2
    return ( (x1 - x2)**2 + (y1 - y2)**2)**0.5

def gradV1x():  # x component of gradV1
    return G*M/r1()**3*x1 + G*m2/r1r2()**3*(x1-x2)
def gradV1y():  # y component of gradV1
    return G*M*m1/r1()**3*y1 + G*m1*m2/r1r2()**3*(y1-y2)

def gradV2x(): # x component of gradV2
    return G*M*m2/r2()**3*x2 + G*m1*m2/r1r2()**3*(x2-x1)
def gradV2y():  # y component of gradV2
    return G*M*m2/r2()**3*y2 + G*m1*m2/r1r2()**3*(y2-y1)

E01 = 0.5*m1*(vx1**2 + vy1**2) - G*M*m1/r1() - G*m1*m2/r1r2() # Total starting energy particle 1
E02 = 0.5*m2*(vx2**2 + vy2**2) - G*M*m2/r2() - G*m1*m2/r1r2() # Total starting energy particle 2

L01 = m1*(x1*vy1 - y1*vx1)   # total angular momentum in start of particle 1, in z-direction
L02 = m2*(x2*vy2 - y2*vx2)   # total angular momentum in start of particle 2, in z-direction

while (time < endtime):

        time = time + dt
        i = i+1

        fx1 = -gradV1x()     # calculating force/mass of particle 1 at time t
        fy1 = -gradV1y()
        vxm1 = vx1 + 0.5*dt*fx1 # first part og verlet algorthm for particle 1
        vym1 = vy1 + 0.5*dt*fy1

        fx2 = -gradV2x()     # calculating force/mass of particle 2 at time t
        fy2 = -gradV2y()
        vxm2 = vx2 + 0.5*dt*fx2 # first part og verlet algorthm for particle 2
        vym2 = vy2 + 0.5*dt*fy2

        x1 = x1 + vxm1*dt       # updating new position, particle 1
        y1 = y1 + vym1*dt
        x2 = x2 + vxm2*dt       # updating new position, particle 2
        y2 = y2 + vym2*dt
        
        fx1 = -gradV1x()     # calculating new force/mass to update new speed, particle 1
        fy1 = -gradV1y()
        vx1 = vxm1 + 0.5*dt*fx1 # calculating new speed, particle 1
        vy1 = vym1 + 0.5*dt*fy1
        
        fx2 = -gradV2x()     # calculating new force/mass to update new speed, particle 2
        fy2 = -gradV2y()
        vx2 = vxm2 + 0.5*dt*fx2 # calculating new speed, particle 2
        vy2 = vym2 + 0.5*dt*fy2

        if (i % plotdensity == 0):  # writing to files only each <plotdensity>'th iteration
            E1 = 0.5*m1*(vx1**2 + vy1**2) - G*M*m1/r1() - G*m1*m2/r1r2() # total energy of particle 1
            E2 = 0.5*m2*(vx2**2 + vy2**2) - G*M*m2/r2() - G*m1*m2/r1r2() # total energy of particle 2
            
            L1 = m1*(x1*vy1 - y1*vx1)   # total angular momentum of particle 1, in z-direction
            L2 = m2*(x2*vy2 - y2*vx2)   # total angular momentum of particle 2, in z-direction
            
            f.write("%f %f\n" % (x1,y1))    # writing to files
            f2.write("%f %f\n" % (x2,y2))
            f3.write("%f %f\n" % (time,E1))
            f4.write("%f %f\n" % (time, E2))
            f5.write("%f %f\n" % (time, L1))
            f6.write("%f %f\n" % (time, L2))

            if E01 > 0 and (abs(1-E1/E01) > maxEnergyDeviation):
                    maxEnergyDeviation = 1 - E1/E01
            if E02 > 0 and (abs(1-E2/E02) > maxEnergyDeviation):
                    maxEnergyDeviation = 1 - E2/E02
            if L01 > 0 and (abs(1-L1/L01) > maxAngularMomentumDeviation):
                    maxAngularMomentumDeviation = L1/L01 - 1
            if L02 > 0 and (abs(1-L2/L02) > maxAngularMomentumDeviation):
                    maxAngularMomentumDeviation = L2/L02 - 1

print("Maximum deviation of energy for a particle relative to start energy is %f\n" % maxEnergyDeviation)
print("Maximum deviation of angular momentum for a particle relative to start angular momentum is %f\n" % maxAngularMomentumDeviation)

#closing files
f.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()

        