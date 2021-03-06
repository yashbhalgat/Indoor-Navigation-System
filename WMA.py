from __future__ import division
from math import sqrt
def lowpass(am, a):   # a is the weight
    n = len(am)
    #mag_field_norm = range(n-1)
    gravity_norm = range(n-1)
    #mag_field = [[0,0,0] for i in range(n-1)]
    gravity = [[0,0,0] for i in range(n-1)]
    a_highpass = [range(n-1),range(n-1),range(n-1)]
    a_highpass_norm = range(n-1)
    gravity[0][0] = am[0][0]
    gravity[0][1] = am[0][1]
    gravity[0][2] = am[0][2]
    #mag_field[0][0] = am[0][3]
    #mag_field[0][1] = am[0][4]
    #mag_field[0][2] = am[0][5]
    gravity_norm[0] = sqrt(gravity[0][0]**2+gravity[0][1]**2+gravity[0][2]**2)
    #mag_field_norm[0] = sqrt(mag_field[0][0]**2+mag_field[0][1]**2+mag_field[0][2]**2)
    for i in range(1,n-1):
        gravity[i][0] = a*gravity[i-1][0]+(1-a)*am[i][0]
        gravity[i][1] = a*gravity[i-1][1]+(1-a)*am[i][1]
        gravity[i][2] = a*gravity[i-1][2]+(1-a)*am[i][2]
        #mag_field[i][0] = a*mag_field[i-1][0]+(1-a)*am[i][3]
        #mag_field[i][1] = a*mag_field[i-1][1]+(1-a)*am[i][4]
        #mag_field[i][2] = a*mag_field[i-1][2]+(1-a)*am[i][5]
        gravity_norm[i] = sqrt(gravity[i][0]**2+gravity[i][1]**2+gravity[i][2]**2)
        #mag_field_norm[i] = sqrt(mag_field[i][0]**2+mag_field[i][1]**2+mag_field[i][2]**2)
    for i in range(1,n-1):
        a_highpass[0][i] = am[i][0]-gravity[i][0]
        a_highpass[1][i] = am[i][1]-gravity[i][1]
        a_highpass[2][i] = am[i][2]-gravity[i][2]
        a_highpass_norm[i] = sqrt(a_highpass[0][i]**2+a_highpass[1][i]**2+a_highpass[2][i]**2)
    cutoff_f = 1000/((a/(1-a))*20)
    print cutoff_f
    return gravity, gravity_norm, a_highpass, a_highpass_norm, n


myfile = open('FFTdata4.txt','r')
am = myfile.read().split('\n')
myfile.close()
am_norm = range(len(am)-1)
#am_norm1 = range(len(am)-1)
for i in range(len(am)-1):
    am[i] = am[i].split(',')
    am[i][0] = float(am[i][0])
    am[i][1] = float(am[i][1])
    am[i][2] = float(am[i][2])
    #am[i][3] = float(am[i][3])
    #am[i][4] = float(am[i][4])
    #am[i][5] = float(am[i][5])
    am_norm[i] = sqrt(am[i][0]**2+am[i][1]**2+am[i][2]**2)
    #am_norm1[i] = sqrt(am[i][3]**2+am[i][4]**2+am[i][5]**2)
gravity, gravity_norm, a_highpass, a_highpass_norm, n = lowpass(am, 0.91)
gravity_transpose = [range(n-1),range(n-1),range(n-1)]
#mag_field_transpose = [range(n-1),range(n-1),range(n-1)]
am_t = [range(n-1),range(n-1),range(n-1)]#,range(n-1),range(n-1),range(n-1)]
for i in range(n-1):
    gravity_transpose[0][i] = gravity[i][0]
    gravity_transpose[1][i] = gravity[i][1]
    gravity_transpose[2][i] = gravity[i][2]
    #mag_field_transpose[0][i] = mag_field[i][0]
    #mag_field_transpose[1][i] = mag_field[i][1]
    #mag_field_transpose[2][i] = mag_field[i][2]
    am_t[0][i] = am[i][0]
    am_t[1][i] = am[i][1]
    am_t[2][i] = am[i][2]
    #am_t[3][i] = am[i][3]
    #am_t[4][i] = am[i][4]
    #am_t[5][i] = am[i][5]
#print gravity_norm[n-2]

def WMA(values, sampling_time, time):
    n = int(time/sampling_time)
    l = len(values)
    values_WMA = range(l)
    for i in range(n):
        values_WMA[i] = values[i]
    for i in range(n,l):
        values_WMA[i] = sum([x for x in [values[i-k] for k in range(n)]])/sum([1 for k in range(n)])
    return values_WMA
gravity_WMA = [[0,0,0] for i in range(n-1)]
gravity_norm_WMA = range(n-1)
gravity_WMA[0] = WMA(gravity_transpose[0], 0.02, 1)
gravity_WMA[1] = WMA(gravity_transpose[1], 0.02, 1)
gravity_WMA[2] = WMA(gravity_transpose[2], 0.02, 1)
am_norm_WMA = WMA(am_norm, 0.033, 2)

import matplotlib.pyplot as plt
plt.subplot(221)
plt.plot(gravity_WMA[0], 'r--', gravity_transpose[0] , 'bs', label = 'abcd')
plt.subplot(222)
plt.plot(gravity_WMA[1], 'r--', gravity_transpose[1] , 'bs')
plt.subplot(223)
plt.plot(gravity_WMA[2], 'r--', gravity_transpose[2] , 'bs')
plt.subplot(224)
plt.plot(am_norm_WMA, 'r--', am_norm , 'bs')
plt.show()
