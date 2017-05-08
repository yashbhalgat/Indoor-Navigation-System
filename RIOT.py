#am is [[ax, ay, az, mx, my, mz] for different time instants] in body frame...... i.e am is a nx6 "matrix"
#am_t is 6xn matrix, i.e it has 6 lists: ax, ay, az, mx, my and mz in body frame each having n components
#am_earth and am_earth_t is in earth frame

from __future__ import division
from numpy import array, sqrt, cross, dot, concatenate
from Euler_from_Quaternion01 import Quat
from numpy.linalg import norm, inv
from math import acos, atan

def QuEST01(v_mag_earth,v_grav_earth,v_mag_mobile,v_grav_mobile, weight):
    # q_BO_est = Estimated quaternion
    # b1,b2 = the 2 measurement vectors in body frames
    # r1,r2 = the 2 measurement vectors in reference frames (orbit/inertial orbit used)
    # b1 is the measurement with very high accurace as compared to b2
    # a1 and a2 are Relative weights of b1 and b2

    #### start TRIAD ####
    r1 = v_mag_earth/norm(v_mag_earth)
    r2 = v_grav_earth/norm(v_grav_earth)
    b1 = v_mag_mobile/norm(v_mag_mobile)
    b2 = v_grav_mobile/norm(v_grav_mobile)

    a1 = weight
    a2 = 1-a1

    #b3
    b3 = cross(b1,b2)
    b3 = b3/norm(b3) 
    #r3
    r3 = cross(r1,r2)
    r3 = r3/norm(r3) 

    #mu
    mu = ( 1 + dot(b3,r3) ) * ( a1*dot(b1,r1) + a2*dot(b2,r2) ) 
    mu =  mu + dot( cross(b3,r3) , ( a1*cross(b1,r1) + a2*cross(b2,r2) ) ) 

    #nu
    nu = dot( ( b3 + r3 ) , ( a1*cross(b1,r1) + a2*cross(b2,r2) )  ) 

    #rho
    rho = sqrt( mu*mu + nu*nu )

    #q
    if( mu > 0):
        k_q_triad = (1/( 2*sqrt( rho*(rho+mu)*( 1 + dot(b3,r3) ) ) ) )
        v_q_triad = ( rho + mu )*( cross(b3,r3) ) + nu*(b3 + r3) 
        s_q_triad = ( rho + mu )* ( 1 + dot(b3,r3) ) 
        q_triad = concatenate((v_q_triad,array([s_q_triad])),axis=0)
        q_triad = k_q_triad * q_triad 
    else:
        k_q_triad = (1/( 2*sqrt( rho*(rho-mu)*( 1 + dot(b3,r3) ) ) ) )
        v_q_triad = ( nu )*( cross(b3,r3) ) + ( rho - mu )*(b3 + r3) 
        s_q_triad = ( nu )* ( 1 + dot(b3,r3) ) 
        q_triad = concatenate((v_q_triad,array([s_q_triad])),axis=0)
        q_triad = k_q_triad * q_triad     


    q_triad = q_triad / sqrt(dot(q_triad,q_triad)) 
    q_BO_est = q_triad
    return q_BO_est

def lowpass(am, cutoff_f):   # a is the weight
    a = (50/cutoff_f)/(1+(50/cutoff_f))
    n = len(am)
    mag_field_norm = range(n-1)
    gravity_norm = range(n-1)
    mag_field = [[0,0,0] for i in range(n-1)]
    gravity = [[0,0,0] for i in range(n-1)]
    omega = [[0,0,0] for i in range(n-1)]
    a_highpass = [range(n-1),range(n-1),range(n-1)]
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
        gravity_norm[i] = sqrt(gravity[i][0]**2+gravity[i][1]**2+gravity[i][2]**2)
        mag_field_norm[i] = sqrt(mag_field[i][0]**2+mag_field[i][1]**2+mag_field[i][2]**2)
    for i in range(1,n-1):
        a_highpass[0][i] = am[i][0]-gravity[i][0]
        a_highpass[1][i] = am[i][1]-gravity[i][1]
        a_highpass[2][i] = am[i][2]-gravity[i][2]
        a_highpass_norm[i] = sqrt(a_highpass[0][i]**2+a_highpass[1][i]**2+a_highpass[2][i]**2)
    #cutoff_f = 1000/((a/(1-a))*20)
    
    print a
    return gravity, gravity_norm, mag_field, mag_field_norm, a_highpass, a_highpass_norm, omega, n

myfile = open('accmaggyro27.txt','r')
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
gravity, gravity_norm, mag_field, mag_field_norm, a_highpass, a_highpass_norm, omega, n = lowpass(am, 5)
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


def RIOT(window, sampling_time):       ### RIOT algo from "Robust pedestrian dead reckoning(R-PDR) for arbitrary mobile device placement"
    no = int(window/sampling_time)
    Q = [[0,0,0,0] for i in range(n-1)]
    gravity_final = [[0,0,0] for i in range(n-1)]
    for i in range(no):
        gravity_final[i] = gravity[i]
        Q[i] = Quat(QuEST01([4.4,38.8,-17.755708694458008],[0,0,1],mag_field[i],gravity_final[i],0.5))
    for i in range(no,n-1):
        s = array([0,0,0])
        for k in range(no):
            s = s+inv(Q[i-1-k]._quat2transform()).dot(array(gravity[i-1-k]).T)
        print "sum ",i," done"
        gravity_final[i] = Q[i-1]._quat2transform().dot((s/no).T)
        print "gravity ",i," done"
#	temp=Q[i-1-k]._quat2transform().dot(array(gravity[i-1-k]).T
#        print temp 
#        print gravity_final[i]
        Q[i] = Quat(QuEST01([41.5, -1.3, -16],[0,0,1],mag_field[i],list(gravity_final[i]),0.5))
    return gravity_final            
sampling_time=0.025
gravity_final = RIOT(2, sampling_time)
print "done"
gravity_final_t = array(gravity_final).T
gravity_final_norm = [norm(gravity_final[i]) for i in range(n-1)]

import matplotlib.pyplot as plt
plt.subplot(221)
plt.plot(gravity_final_t[0], 'r--', gravity_transpose[0], 'b--')
plt.subplot(222)
plt.plot(gravity_final_t[1], 'r--', gravity_transpose[1], 'b--')
plt.subplot(223)
plt.plot(gravity_final_t[2], 'r--', gravity_transpose[2], 'b--')
plt.subplot(224)
plt.plot(gravity_final_norm, 'r--', gravity_norm, 'b--')
plt.show()

