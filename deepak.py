#am is [[ax, ay, az, mx, my, mz, gyrox, gyroy, gyroz, time] for different time instants] in body frame...... i.e am is a nx6 "matrix"
#am_t is 6xn matrix, i.e it has 10 lists: ax, ay, az, mx, my, mz, gyrox, gyroy, gyroz, time in body frame each having n components
#am_earth and am_earth_t are the accelerations in earth frame

from __future__ import division
from numpy import *
from Euler_from_Quaternion01 import Quat
from numpy.linalg import norm, inv
from math import acos, atan2
from PCA import princomp
from Quaternion import QuEST01
from scipy import hamming

def stft(x, fs, framesz):
    framesamp = int(framesz*fs)
    hopsamp = 1
    w = hamming(framesamp)
    X = array([fft.fft(w*x[i:i+framesamp]) 
                     for i in range(0, len(x)-framesamp, hopsamp)])
    return X

def lowpass(x, cutoff_f): # x is a list of "n" 3 element lists
    a = (50/cutoff_f)/(1+(50/cutoff_f))
    n = len(x)
    x_lp = [[0,0,0] for i in range(n-1)]
    x_lp[0] = x[0]
    x_lp_norm = range(n-1)
    x_lp_norm[0] = sqrt(x_lp[0][0]**2+x_lp[0][1]**2+x_lp[0][2]**2)
    for i in range(1,n-1):
        x_lp[i][0] = a*x_lp[i-1][0]+(1-a)*x[i][0]
        x_lp[i][1] = a*x_lp[i-1][1]+(1-a)*x[i][1]
        x_lp[i][2] = a*x_lp[i-1][2]+(1-a)*x[i][2]
        x_lp_norm[i] = sqrt(x_lp[i][0]**2+x_lp[i][1]**2+x_lp[i][2]**2)
    return x_lp, x_lp_norm

def matrix_correct(r):
    r1 = array(empty((3, 3), float))
    error = r[0].dot(r[1].T)
    r1[0] = r[0] - (error/2.0)*r[1]
    r1[1] = r[1] - (error/2.0)*r[0]
    r1[2] = cross(r1[0],r1[1])
    r1[0] = r1[0]/norm(r1[0])
    r1[1] = r1[1]/norm(r1[1])
    r1[2] = r1[2]/norm(r1[2])
    return r1

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
    k = abs(c[i:j]).argmax()+i
    return f[k]

def Quat_CF(gamma, OMEGA_q, Quaternion, am, Time):     ###ref: eqn [33] and [41] from paper "An_efficient_orientation_filter_for_inertial_and_inertialmagnetic_sensor_arrays.pdf"
    n = len(am)
    Q_CF = [[0,0,0,1] for i in range(n-1)]
    Q_CF[0] = Quaternion[0].q
    for i in range(1,n-1):
        Q_CF[i] = gamma*(Quaternion[i].q-Q_CF[i-1]) + (1-gamma)*(Q_CF[i-1]+0.5*array(OMEGA_q).dot(array(Q_CF[i-1]).T)*(Time[i]-Time[i-1]))
    return Q_CF

def PCA_dynamic(window, sampling_time, a_earth_t, am):
    n = len(am)
    no = int(window/sampling_time)
    heading_b = [range(no, n-1), range(no, n-1), range(no, n-1)]
    heading_angle = range(no, n-1)
    for i in range(n-no-2):
        A = array([a_earth_t[0][i:no+i],a_earth_t[1][i:no+i]])#,a_earth_t[2][i:no+i]])
        coeff, score, latent = princomp(A.T)
        if(sum(A[0])/len(A[0])>0):
        	heading_b[0][i] = abs(coeff[0][latent.argmax()])
        elif(sum(A[0])/len(A[0])<=0):
        	heading_b[0][i] = -abs(coeff[0][latent.argmax()])
        if(sum(A[1])/len(A[1])>0):
        	heading_b[1][i] = abs(coeff[1][latent.argmax()])
        elif(sum(A[1])/len(A[1])<=0):
        	heading_b[1][i] = -abs(coeff[1][latent.argmax()])
        heading_angle[i] = atan2(heading_b[1][i],heading_b[0][i])
        #heading_b[2][i] = abs(coeff[2][latent.argmax()])
        #print coeff[0][latent.argmax()], coeff[1][latent.argmax()], coeff[2][latent.argmax()]
        #print coeff[0][latent.argmax()], coeff[1][latent.argmax()], latent
    return heading_b, heading_angle

def main():
	myfile = open('final3.txt','r')
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
		am_norm[i] = sqrt(am[i][0]**2+am[i][1]**2+am[i][2]**2)
		am_norm1[i] = sqrt(am[i][3]**2+am[i][4]**2+am[i][5]**2)

	gravity, gravity_norm = lowpass([x[0:3] for x in am], 1)
	mag_field, mag_field_norm = lowpass([x[3:6] for x in am], 1)
	omega = [[0,0,0] for i in range(n-1)]
	for i in range(n-1):    
		omega[i][0] = am[i][6]
		omega[i][1] = am[i][7]
		omega[i][2] = am[i][8]

	gravity_transpose = array(gravity).T
	mag_field_transpose = array(mag_field).T
	am_t = array(am[:-1]).T

	#roll, pitch, yaw, Quaternion and all... just in case, you want to use it
	roll, pitch, yaw = range(n-1),range(n-1),range(n-1)
	Quaternion = [[0,0,0,0] for i in range(n-1)]
	theta = range(n-1)
	e = [range(n-1),range(n-1),range(n-1)]
	for i in range(n-1):
		Quaternion[i] = Quat(QuEST01([22.9, -55.8, -17.8],[0,0,1],mag_field[i],gravity[i],0.5))
		roll[i], pitch[i], yaw[i] = Quaternion[i].roll, Quaternion[i].dec, Quaternion[i].ra 
		#print roll[i], pitch[i], yaw[i]
		if(acos(Quaternion[i].q[3])==0):
		    theta[i] = 0.01
		else: theta[i] = 2.0*acos(Quaternion[i].q[3])
		e[0][i] = Quaternion[i].q[0]/sin(theta[i]/2)
		e[1][i] = Quaternion[i].q[1]/sin(theta[i]/2)
		e[2][i] = Quaternion[i].q[2]/sin(theta[i]/2)

	rmat = [empty((3, 3), float) for i in range(n-1)]

	Time = am_t[9]
	Time=append(Time, am_t[9][-2])
	dT = range(n-1)
	for i in range(n-1):
		dT[i] = Time[i+1]-Time[i]
	dT1 = []
	for x in dT:
		if (x!=0):
		    dT1.append(x)

	OMEGA_q = [empty((4, 4), float) for i in range(n-1)]     ### the 4x4 matrix related to the rate of change of quaternion
	for i in range(n-1):
		OMEGA_q = [[0, -omega[i][2], omega[i][1], omega[i][0]],[omega[i][2], 0, -omega[i][0], omega[i][1]],[-omega[i][1], omega[i][0], 0, omega[i][2]],[-omega[i][0], -omega[i][1], -omega[i][2], 0]]

	rmat_gyro = [empty((3, 3), float) for i in range(n-1)]

	Q_CF = Quat_CF(0.01, OMEGA_q, Quaternion, am, Time)

	a_earth = [[0,0,0] for i in range(n-1)]
	gravity_earth = [[0,0,0] for i in range(n-1)]
	for i in range(n-1):
		a_earth[i] = (Quat(Q_CF[i]/norm(Q_CF[i]))._quat2transform()).dot(array(am[i][0:3]).T)
		gravity_earth[i] = (Quat(Q_CF[i]/norm(Q_CF[i]))._quat2transform()).dot(array(gravity[i]).T)
		#print a_earth[i]

	a_earth_t = array(a_earth).T
	
	sampling_time=0.025
	fs, framesz = 1//sampling_time, 2
	no = int(framesz/sampling_time)
	
	heading_b, heading_angle = PCA_dynamic(framesz, sampling_time, a_earth_t, am)
	
	a = am_norm
	sp = stft(a, fs, framesz)
	freq = fft.fftfreq(int(fs*framesz),d=sampling_time)

	max_f = []
	for x in sp:
		max_f.append(max_freq(x, freq, 0.8, 3))

	'''import matplotlib.pyplot as plt
	plt.subplot(221)
	plt.plot(sampling_time*linspace(1, len(max_f)-1, len(max_f)), max_f, '-ro')
	plt.subplot(222)
	plt.plot(a, '-ro')
	plt.show()'''
	
	h = 1.75
	step_length = [(0.15+0.14*f)*h for f in max_f]
	x=[0]
	y=[0]
	for i in range(len(max_f[:-2])):
		x.append(x[-1]+(0.15+0.14*f)*h*f*cos(heading_angle[i])*dT[no:][i])
		y.append(y[-1]+(0.15+0.14*f)*h*f*sin(heading_angle[i])*dT[no:][i])
	
	return [[i+1, dT[no:][i], step_length[:-2][i], max_f[:-2][i], array(heading_b).T[i]] for i in range(n-no-1)]


