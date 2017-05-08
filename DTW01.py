import mlpy

import matplotlib.cm as cm
from numpy import *

x = [sin(2*pi*t) for t in linspace(0,5,100)]
y = [sin(2*pi*t)+random.random() for t in linspace(0,5,100)]
dist, cost, path = mlpy.dtw_std(x, y, dist_only=False)
#dist
#0.0
print dist
NewX=[]
NewY=[]
for i in range(len(path[0])):
	NewX.append(x[path[0][i]])
for j in range(len(path[1])):
	NewY.append(y[path[1][j]])
dist1, cost1, path1 = mlpy.dtw_std(NewX, NewY, dist_only=False)
print dist1
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt1
fig = plt.figure(1)
ax = fig.add_subplot(111)
plot1 = plt.imshow(cost.T, origin='lower', cmap=cm.gray, interpolation='nearest')
plot2 = plt.plot(path[0], path[1], 'w')
xlim = ax.set_xlim((-0.5, cost.shape[0]-0.5))
ylim = ax.set_ylim((-0.5, cost.shape[1]-0.5))
plt.show()
plt1.subplot(221)
plt1.plot(x,'r--')
plt1.subplot(222)
plt1.plot(y,'r--')
plt1.subplot(223)
plt1.plot(NewX,'r--')
plt1.subplot(224)
plt1.plot(NewY,'r--')
plt1.show()
