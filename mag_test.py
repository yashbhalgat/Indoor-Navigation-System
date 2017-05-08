from numpy import array
myfile = open('mag_test2.txt','r')
am = myfile.read().split('\n')
myfile.close()
am_norm = range(len(am)-1)
am_norm1 = range(len(am)-1)
n = len(am)
for i in range(len(am)-1):
    am[i] = am[i].split(',')
    am[i][0] = float(am[i][0])
    am[i][1] = float(am[i][1])
    am[i][2] = float(am[i][2])
    am[i][3] = float(am[i][3])
    am[i][4] = float(am[i][4])
    am[i][5] = float(am[i][5])
    am[i][6] = float(am[i][6])
    am[i][7] = float(am[i][7])
    am[i][8] = float(am[i][8])
    am[i][9] = float(am[i][9])    
    
am_t = array(am[:-1]).T
mag = [am_t[3], am_t[4], am_t[5]]
print mag
from matplotlib import pyplot
import pylab
from mpl_toolkits.mplot3d import Axes3D
import random

fig = pylab.figure()
ax = Axes3D(fig)
ax.set_xlim([-100,100])
ax.set_ylim([-100,100])
ax.set_zlim([-100,100])
ax.scatter(mag[0], mag[1], mag[2])
pylab.show()

