import stft
import numpy
T = 0.05
fs = 5000.0
nsamples = T * fs
t = numpy.linspace(0, T, nsamples, endpoint=False)
a = 0.02
f0 = 600.0
x = 0.1 * numpy.sin(2 * numpy.pi * 1.2 * numpy.sqrt(t))
x += 0.01 * numpy.cos(2 * numpy.pi * 312 * t + 0.1)
x += a * numpy.cos(2 * numpy.pi * f0 * t + .11)
x += 0.03 * numpy.cos(2 * numpy.pi * 2000 * t)
import matplotlib.pyplot as plt
plt.plot(x, 'r--')
plt.show()
print stft.process(x, scipy.signal.get_window('hamming',100))
