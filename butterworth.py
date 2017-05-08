from scipy.signal import butter, lfilter
from math import sqrt
import matplotlib.pyplot as plt
import numpy as np

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y
'''
def lowpass(am, a):   # a is the weight
    n = len(am)
    mag_field_norm = range(n-1)
    gravity_norm = range(n-1)
    mag_field = [[0,0,0] for i in range(n-1)]
    gravity = [[0,0,0] for i in range(n-1)]
    a_highpass = [range(n-1),range(n-1),range(n-1)]
    a_highpass_norm = range(n-1)
    gravity[0][0] = am[0][0]
    gravity[0][1] = am[0][1]
    gravity[0][2] = am[0][2]
    mag_field[0][0] = am[0][3]
    mag_field[0][1] = am[0][4]
    mag_field[0][2] = am[0][5]
    gravity_norm[0] = sqrt(gravity[0][0]**2+gravity[0][1]**2+gravity[0][2]**2)
    mag_field_norm[0] = sqrt(mag_field[0][0]**2+mag_field[0][1]**2+mag_field[0][2]**2)
    for i in range(1,n-1):
        gravity[i][0] = a*gravity[i-1][0]+(1-a)*am[i][0]
        gravity[i][1] = a*gravity[i-1][1]+(1-a)*am[i][1]
        gravity[i][2] = a*gravity[i-1][2]+(1-a)*am[i][2]
        mag_field[i][0] = a*mag_field[i-1][0]+(1-a)*am[i][3]
        mag_field[i][1] = a*mag_field[i-1][1]+(1-a)*am[i][4]
        mag_field[i][2] = a*mag_field[i-1][2]+(1-a)*am[i][5]
        gravity_norm[i] = sqrt(gravity[i][0]**2+gravity[i][1]**2+gravity[i][2]**2)
        mag_field_norm[i] = sqrt(mag_field[i][0]**2+mag_field[i][1]**2+mag_field[i][2]**2)
    for i in range(1,n-1):
        a_highpass[0][i] = am[i][0]-gravity[i][0]
        a_highpass[1][i] = am[i][1]-gravity[i][1]
        a_highpass[2][i] = am[i][2]-gravity[i][2]
        a_highpass_norm[i] = sqrt(a_highpass[0][i]**2+a_highpass[1][i]**2+a_highpass[2][i]**2)
    cutoff_f = 1000/((a/(1-a))*40)
    print cutoff_f
    return gravity, gravity_norm, mag_field, mag_field_norm, a_highpass, a_highpass_norm, n
'''

myfile = open('walking_mag_static23.txt','r')
am = myfile.read().split('\n')
myfile.close()
n = len(am)
am_norm = range(len(am)-1)
am_norm1 = range(len(am)-1)
for i in range(len(am)-1):
    am[i] = am[i].split(',')
    am[i][0] = float(am[i][0])
    am[i][1] = float(am[i][1])
    am[i][2] = float(am[i][2])
    am[i][3] = float(am[i][3])
    am[i][4] = float(am[i][4])
    am[i][5] = float(am[i][5])
    am_norm[i] = sqrt(am[i][0]**2+am[i][1]**2+am[i][2]**2)
    am_norm1[i] = sqrt(am[i][3]**2+am[i][4]**2+am[i][5]**2)

#gravity_transpose = [range(n-1),range(n-1),range(n-1)]

am_t = [range(n-1),range(n-1),range(n-1),range(n-1),range(n-1),range(n-1)]
for i in range(n-1):
    am_t[0][i] = am[i][0]
    am_t[1][i] = am[i][1]
    am_t[2][i] = am[i][2]
    am_t[3][i] = am[i][3]
    am_t[4][i] = am[i][4]
    am_t[5][i] = am[i][5]

gravity = [0,0,0]
fs = 17
gravity[0] = butter_bandpass_filter(am_t[0], 0, 5, fs, order=3)
gravity[1] = butter_bandpass_filter(am_t[1], 0, 5, fs, order=3)
gravity[2] = butter_bandpass_filter(am_t[2], 0, 5, fs, order=3)
gravity_norm = butter_bandpass_filter(am_norm, 0, 5, fs, order=3)

a = range(n-1)
a = gravity_norm
sp = np.fft.fft(a)
freq = np.fft.fftfreq(len(a),d=1.0/fs)

i=0
j=0

for x in freq:
    if(abs(x-0.5)<=0.2): i=list(freq).index(x)
    if(abs(x-3)<=0.2): j=list(freq).index(x)

k = abs(sp[i:j]).argmax()+i
print freq[k]



plt.subplot(321)
plt.plot(am_t[0], 'r--', gravity[0] , 'b--')
plt.subplot(322)
plt.plot(am_t[1], 'r--', gravity[1] , 'b--')
plt.subplot(323)
plt.plot(am_t[2], 'r--', gravity[2] , 'b--')
plt.subplot(324)
plt.plot(am_norm, 'r--', gravity_norm , 'b--')
#plt.plot(freq, sp.real,'r--', freq, sp.imag, 'r--')
plt.subplot(325)
plt.plot(freq, abs(sp), '-ro')
#plt.plot(am_t[1], 'r--')
plt.show()

plt.show()
