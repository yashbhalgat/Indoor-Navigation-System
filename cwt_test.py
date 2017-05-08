import numpy as np
from scipy import signal
sig = [5*np.sin(2*5*3.142*t)+3*np.sin(2*10*3.142*t) for t in np.linspace(0,2,100)]
wavelet = signal.ricker
widths = np.arange(1, 11)
cwtmatr = signal.cwt(sig, wavelet, widths)
