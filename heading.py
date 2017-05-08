#am is [[ax, ay, az, mx, my, mz] for different time instants] in body frame...... i.e am is a nx6 "matrix"
#am_t is 6xn matrix, i.e it has 6 lists: ax, ay, az, mx, my and mz in body frame each having n components
#am_earth and am_earth_t is in earth frame



from __future__ import division
from numpy import *
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

def princomp(A):
 """ performs principal components analysis 
     (PCA) on the n-by-p data matrix A
     Rows of A correspond to observations, columns to variables. 

 Returns :  
  coeff :
    is a p-by-p matrix, each column containing coefficients 
    for one principal component.
  score : 
    the principal component scores; that is, the representation 
    of A in the principal component space. Rows of SCORE 
    correspond to observations, columns to components.
  latent : 
    a vector containing the eigenvalues 
    of the covariance matrix of A.
 """
 # computing eigenvalues and eigenvectors of covariance matrix
 M = (A-mean(A.T,axis=1)).T # subtract the mean (along columns)
 [latent,coeff] = linalg.eig(cov(M)) # attention:not always sorted
 score = dot(coeff.T,M) # projection of the data in the new space
 return coeff,score,latent

def lowpass(am, cutoff_f):   # a is the weight
    a = (50/cutoff_f)/(1+(50/cutoff_f))
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
    #cutoff_f = 1000/((a/(1-a))*20)
    
    print a
    return gravity, gravity_norm, mag_field, mag_field_norm, a_highpass, a_highpass_norm, n

myfile = open('walking_mag_static9.txt','r')
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
    am_norm[i] = sqrt(am[i][0]**2+am[i][1]**2+am[i][2]**2)
    am_norm1[i] = sqrt(am[i][3]**2+am[i][4]**2+am[i][5]**2)
gravity, gravity_norm, mag_field, mag_field_norm, a_highpass, a_highpass_norm, n = lowpass(am, 4)
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

roll, pitch, yaw = range(n-1),range(n-1),range(n-1)
Quaternion = [[0,0,0,0] for i in range(n-1)]
theta = range(n-1)
e = [range(n-1),range(n-1),range(n-1)]
for i in range(n-1):
    Quaternion[i] = Quat(QuEST01([4.4,38.8,-17.755708694458008],[0,0,1],mag_field[i],gravity[i],0.5))
    roll[i], pitch[i], yaw[i] = Quaternion[i].roll, Quaternion[i].dec, Quaternion[i].ra 
    #print roll[i], pitch[i], yaw[i]
    if(acos(Quaternion[i].q[3])==0):
        theta[i] = 0.01
    else: theta[i] = 2.0*acos(Quaternion[i].q[3])
    e[0][i] = Quaternion[i].q[0]/sin(theta[i]/2)
    e[1][i] = Quaternion[i].q[1]/sin(theta[i]/2)
    e[2][i] = Quaternion[i].q[2]/sin(theta[i]/2)

a_earth = [[0,0,0] for i in range(n-1)]
for i in range(n-1):
    a_earth[i] = (Quaternion[i]._quat2transform()).dot(array(am[i][0:3]).T)
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

import matplotlib.pyplot as plt

plt.subplot(221)
plt.plot(heading_b[0], 'r--')
plt.title('heading_b_x')
plt.ylim(ymax = 2, ymin = -2)
plt.subplot(222)
plt.plot(heading_b[1], 'r--')
plt.title('heading_b_y')
plt.ylim(ymax = 2, ymin = -2)
'''plt.subplot(223)
plt.plot(heading_b[2], 'r--')
plt.title('heading_b_z')
plt.ylim(ymax = 2, ymin = -2)'''
plt.subplot(223)
plt.plot(heading_angle, 'r--')
plt.title('heading_angle')
#plt.ylim(ymax = 2, ymin = -2)
plt.show()
