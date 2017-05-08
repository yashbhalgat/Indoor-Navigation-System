import numpy as np
import matplotlib.pyplot as plt
import math as m
def msin(i):
    return np.sin(i*np.pi*10)
a = map(msin,[0.01*i for i in range(500)])#[np.sin(i*pi*10) for i in np.linspace(0,5,100)]
#plt.plot(a)
sp = np.fft.fft(a)
freq = np.fft.fftfreq(len(a),d=0.01)

plt.plot(freq, sp.real, freq, sp.imag)

plt.show()
