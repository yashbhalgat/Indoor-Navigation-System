import matplotlib.pyplot as plt
from cmath import exp, pi
from math import *
import numpy as np
import scipy

def stft(x, fs, framesz, hop):
    framesamp = int(framesz*fs)
    hopsamp = int(hop*fs)
    w = scipy.hamming(framesamp)
    X = scipy.array([np.fft.fft(w*x[i:i+framesamp]) 
                     for i in range(0, len(x)-framesamp, hopsamp)])
    return X

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
    #print cutoff_f
    return gravity, gravity_norm, a_highpass, a_highpass_norm, n


myfile = open('STFT_testing1.txt','r')
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
gravity, gravity_norm, a_highpass, a_highpass_norm, n = lowpass(am, 0.99)
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

def WMA(values, sampling_time, time):
    n = int(time/sampling_time)
    l = len(values)
    values_WMA = range(l)
    for i in range(n):
        values_WMA[i] = values[i]
    for i in range(n,l):
        values_WMA[i] = sum([x for x in [values[i-k] for k in range(n)]])/sum([1 for k in range(n)])
    return values_WMA
a_highpass_normWMA = WMA(a_highpass_norm, 0.33, 2)

sampling_time=0.025
fs, framesz, hop = 1//sampling_time, 6, 0.05
a = range(n-1)
t = np.linspace(0, 2*3.142*20*5, 5*fs, endpoint=False)

a = am_norm
sp = stft(a, fs, framesz, hop)
freq = np.fft.fftfreq(int(fs*framesz),d=sampling_time)

def max_freq(c, f, lowcut, highcut):   #gives freq in "f" between "lowcut" and "highcut" which has max corresponding coeff in "c"
    i=0
    j=0
    for x in f[0:len(f)//2]:
        if(x<=lowcut):
            i=i+1
        #if(abs(x-lowcut)<=0.2): i=list(f).index(x)
        if(x<=highcut):
            j=j+1
        #if(abs(x-highcut)<=0.2): j=list(f).index(x)
    j=j+1
    print i, j
    k = abs(c[i:j]).argmax()+i
    print abs(c[i:j]), f[i:j]
    return f[k]
print freq
max_f = []
for x in sp:
    max_f.append(max_freq(x, freq, 0.8, 3))


#freq = np.fft.fftfreq(len(a),d=0.033)

'''i=0
j=0

for x in freq:
    if(abs(x-0.5)<=0.2): i=list(freq).index(x)
    if(abs(x-3)<=0.2): j=list(freq).index(x)

k = abs(sp[i:j]).argmax()+i
print freq[k]
'''

#plt.plot(freq, sp.real,'r--', freq, sp.imag, 'r--')
plt.subplot(221)
plt.plot(sampling_time*np.linspace(1, len(max_f)-1, len(max_f)), max_f, '-ro')
plt.subplot(222)
plt.plot(a, '-ro')
#plt.plot(am_t[1], 'r--')
plt.show()

'''
if __name__ == '__main__':
    f0 = 40         # Compute the STFT of a 440 Hz sinusoid
    fs = 80       # sampled at 8 kHz
    T = 5            # lasting 5 seconds
    framesz = 0.050  # with a frame size of 50 milliseconds
    hop = 0.020      # and hop size of 20 milliseconds.

    # Create test signal and STFT.
    t = scipy.linspace(0, T, T*fs, endpoint=False)
    x = scipy.sin(2*scipy.pi*f0*t)
    X = stft(x, fs, framesz, hop)

    # Plot the magnitude spectrogram.
    #plt.figure()
    #plt.imshow(scipy.absolute(X.T), origin='lower', aspect='auto',
                # interpolation='nearest')
    #plt.xlabel('Time')
    #plt.ylabel('Frequency')
    plt.plot( x, 'r--')
    plt.show()
    '''
