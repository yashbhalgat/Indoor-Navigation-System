#am is [[ax, ay, az, mx, my, mz] for different time instants] in body frame...... i.e am is a nx6 "matrix"
#am_t is 6xn matrix, i.e it has 6 lists: ax, ay, az, mx, my and mz in body frame each having n components
#am_earth and am_earth_t is in earth frame

from __future__ import division
from numpy import *
from Euler_from_Quaternion01 import Quat
from numpy.linalg import norm, inv
from math import acos, atan
from PCA import princomp
from Quaternion import QuEST01

def lowpass(am, cutoff_f):   # a is the weight
    a = (50/cutoff_f)/(1+(50/cutoff_f))
    n = len(am)
    mag_field_norm = range(n-1)
    gravity_norm = range(n-1)
    mag_field = [[0,0,0] for i in range(n-1)]
    gravity = [[0,0,0] for i in range(n-1)]
    omega = [[0,0,0] for i in range(n-1)]
    omega_lp = [[0,0,0] for i in range(n-1)]
    a_highpass = [range(n-1),range(n-1),range(n-1)]
    omega_hp = [[0,0,0] for i in range(n-1)]
    a_highpass_norm = range(n-1)
    gravity[0][0] = am[0][0]
    gravity[0][1] = am[0][1]
    gravity[0][2] = am[0][2]
    mag_field[0][0] = am[0][3]
    mag_field[0][1] = am[0][4]
    mag_field[0][2] = am[0][5]
    omega[0][0] = am[0][6]
    omega[0][1] = am[0][7]
    omega[0][2] = am[0][8]
    omega_lp[0] = omega[0]
    gravity_norm[0] = sqrt(gravity[0][0]**2+gravity[0][1]**2+gravity[0][2]**2)
    mag_field_norm[0] = sqrt(mag_field[0][0]**2+mag_field[0][1]**2+mag_field[0][2]**2)
    for i in range(1,n-1):
        gravity[i][0] = a*gravity[i-1][0]+(1-a)*am[i][0]
        gravity[i][1] = a*gravity[i-1][1]+(1-a)*am[i][1]
        gravity[i][2] = a*gravity[i-1][2]+(1-a)*am[i][2]
        mag_field[i][0] = a*mag_field[i-1][0]+(1-a)*am[i][3]
        mag_field[i][1] = a*mag_field[i-1][1]+(1-a)*am[i][4]
        mag_field[i][2] = a*mag_field[i-1][2]+(1-a)*am[i][5]
        omega[i][0] = am[i][6]
        omega[i][1] = am[i][7]
        omega[i][2] = am[i][8]
        omega_lp[i][0] = a*omega_lp[i-1][0]+(1-a)*omega[i][0]
        omega_lp[i][1] = a*omega_lp[i-1][1]+(1-a)*omega[i][1]
        omega_lp[i][2] = a*omega_lp[i-1][2]+(1-a)*omega[i][2]
        gravity_norm[i] = sqrt(gravity[i][0]**2+gravity[i][1]**2+gravity[i][2]**2)
        mag_field_norm[i] = sqrt(mag_field[i][0]**2+mag_field[i][1]**2+mag_field[i][2]**2)
    for i in range(1,n-1):
        a_highpass[0][i] = am[i][0]-gravity[i][0]
        a_highpass[1][i] = am[i][1]-gravity[i][1]
        a_highpass[2][i] = am[i][2]-gravity[i][2]
        a_highpass_norm[i] = sqrt(a_highpass[0][i]**2+a_highpass[1][i]**2+a_highpass[2][i]**2)
        omega_hp[i][0] = omega[i][0]-omega_lp[i][0]
        omega_hp[i][1] = omega[i][1]-omega_lp[i][1]
        omega_hp[i][2] = omega[i][2]-omega_lp[i][2]
    #cutoff_f = 1000/((a/(1-a))*20)
    
    print a
    return gravity, gravity_norm, mag_field, mag_field_norm, a_highpass, a_highpass_norm, omega, omega_hp, n

myfile = open('time_test.txt','r')
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
gravity, gravity_norm, mag_field, mag_field_norm, a_highpass, a_highpass_norm, omega, omega_hp, n = lowpass(am, 3)
gravity_transpose = [range(n-1),range(n-1),range(n-1)]
mag_field_transpose = [range(n-1),range(n-1),range(n-1)]
am_t = [range(n-1),range(n-1),range(n-1),range(n-1),range(n-1),range(n-1)]
for i in range(n-1):
    gravity_transpose[0][i] = gravity[i][0]
    gravity_transpose[1][i] = gravity[i][1]
    gravity_transpose[2][i] = gravity[i][2]
    mag_field_transpose[0][i] = mag_field[i][0]
    mag_field_transpose[1][i] = mag_field[i][1]
    mag_field_transpose[2][i] = mag_field[i][2]
    am_t[0][i] = am[i][0]
    am_t[1][i] = am[i][1]
    am_t[2][i] = am[i][2]
    am_t[3][i] = am[i][3]
    am_t[4][i] = am[i][4]
    am_t[5][i] = am[i][5]
    
#roll, pitch, yaw, e and theta for various readings
roll, pitch, yaw = range(n-1),range(n-1),range(n-1)
Quaternion = [[0,0,0,1] for i in range(n-1)]
theta = range(n-1)
e = [range(n-1),range(n-1),range(n-1)]
for i in range(n-1):
    Quaternion[i] = Quat(QuEST01([41.5, -1.3, -16],[0,0,1],mag_field[i],gravity[i],0.5))
    roll[i], pitch[i], yaw[i] = Quaternion[i].roll, Quaternion[i].dec, Quaternion[i].ra 
    #print roll[i], pitch[i], yaw[i]
    if(acos(Quaternion[i].q[3])==0):
        theta[i] = 0.01
    else: theta[i] = 2.0*acos(Quaternion[i].q[3])
    e[0][i] = Quaternion[i].q[0]/sin(theta[i]/2)
    e[1][i] = Quaternion[i].q[1]/sin(theta[i]/2)
    e[2][i] = Quaternion[i].q[2]/sin(theta[i]/2)

rmat = [empty((3, 3), float) for i in range(n-1)]

Time = range(n)
for i in range(n-1):
    Time[i] = (am[i][9])
Time[n-1] = am[n-2][9]
dT = range(n-1)
for i in range(n-1):
    dT[i] = Time[i+1]-Time[i]
dT1 = []
for x in dT:
    if (x!=0):
        dT1.append(x)

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

OMEGA_q = [empty((4, 4), float) for i in range(n-1)]     ### the 4x4 matrix related to the rate of change of quaternion
for i in range(n-1):
    OMEGA_q = [[0, -omega[i][2], omega[i][1], omega[i][0]],[omega[i][2], 0, -omega[i][0], omega[i][1]],[-omega[i][1], omega[i][0], 0, omega[i][2]],[-omega[i][0], -omega[i][1], -omega[i][2], 0]]


rmat_gyro = [empty((3, 3), float) for i in range(n-1)]
def using_gyro(sampling_time):
    OMEGA = [empty((3, 3), float) for i in range(n-1)]
    sum = 0
    # rmat_gyro = [empty((3, 3), float) for i in range(n-1)]
    rmat_gyro[0] = rmat[0]
    for i in range(1,n-1):
        OMEGA[i] = array([[0, omega[i][2], -omega[i][1]],[-omega[i][2], 0, omega[i][0]],[omega[i][1], -omega[i][0], 0]])
        rmat_gyro[i] = array(matrix(rmat_gyro[i-1])*(eye(3)-matrix(OMEGA[i])*(Time[i]-Time[i-1])))
        sum = sum + omega_hp[i][2]*sampling_time
        rmat_gyro[i] = matrix_correct(rmat_gyro[i])
    Q = Quat(matrix(rmat_gyro[0])*inv(matrix(rmat_gyro[n-2])))
    print sum*180/3.142
    print 2*acos(Q.q[3])*180/3.142#return rmat_gyro

def Quat_CF(gamma):     ###ref: eqn [33] and [41] from paper "An_efficient_orientation_filter_for_inertial_and_inertialmagnetic_sensor_arrays.pdf"
    Q_CF = [[0,0,0,1] for i in range(n-1)]
    Q_CF[0] = Quaternion[0].q
    for i in range(1,n-1):
        Q_CF[i] = gamma*(Quaternion[i].q-Q_CF[i-1]) + (1-gamma)*(Q_CF[i-1]+0.5*array(OMEGA_q).dot(array(Q_CF[i-1]).T)*(Time[i]-Time[i-1]))
    return Q_CF

Q_CF = Quat_CF(0.01)

a_earth = [[0,0,0] for i in range(n-1)]
for i in range(n-1):
    a_earth[i] = (Quat(Q_CF[i]/norm(Q_CF[i]))._quat2transform()).dot(array(am[i][0:3]).T)
    #print a_earth[i]

a_earth_t = array(a_earth).T

def PCA_dynamic(window, sampling_time):
    no = int(window/sampling_time)
    heading_b = [range(no, n-1), range(no, n-1), range(no, n-1)]
    heading_angle = range(no, n-1)
    for i in range(n-no-2):
        A = array([a_earth_t[0][i:no+i],a_earth_t[1][i:no+i]])#,a_earth_t[2][i:no+i]])
        coeff, score, latent = princomp(A.T)
        heading_b[0][i] = abs(coeff[0][latent.argmax()])
        heading_b[1][i] = abs(coeff[1][latent.argmax()])
        heading_angle[i] = atan(heading_b[0][i]/ heading_b[1][i])*180/pi
        #heading_b[2][i] = abs(coeff[2][latent.argmax()])
        #print coeff[0][latent.argmax()], coeff[1][latent.argmax()], coeff[2][latent.argmax()]
        #print coeff[0][latent.argmax()], coeff[1][latent.argmax()], latent
    return heading_b, heading_angle


heading_b, heading_angle = PCA_dynamic(2, 0.02)

#print 2*acos(Quat(Q_CF[0]/norm(Q_CF[0])).q[3])*180/3.142, 2*acos(Quat(Q_CF[n-3]/norm(Q_CF[n-3])).q[3])*180/3.142
#print 2*acos((Quaternion[0].q/norm(Quaternion[0].q))[3])*180/3.142, 2*acos((Quaternion[n-3].q/norm(Quaternion[n-3].q))[3])*180/3.142
#using_gyro(0.02451)
print rmat_gyro[0]
print rmat_gyro[n-2]
print "mean", sum(dT)/len(dT), "std", std(dT1) 

omega_hp_t = array(omega_hp).T
import matplotlib.pyplot as plt
plt.subplot(321)
plt.plot(dT, 'r--')
plt.subplot(322)
plt.plot(omega_hp_t[1], 'r--')
plt.subplot(323)
plt.plot(omega_hp_t[2], 'r--')
plt.subplot(324)
plt.plot([Quat(Q_CF[i]/norm(Q_CF[i])).ra for i in range(n-1)], 'r--')
plt.subplot(325)
plt.plot([Quat(Q_CF[i]/norm(Q_CF[i])).dec for i in range(n-1)], 'r--')
plt.subplot(326)
plt.plot([Quat(Q_CF[i]/norm(Q_CF[i])).roll for i in range(n-1)], 'r--')
plt.show()

