from __future__ import division

from numpy import concatenate, cross, dot, array
from numpy.linalg import norm
from math import *

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

if __name__=="__main__":
    myfile = open('test_cases_for_Quest.csv','r')
    s = myfile.read().split('\n')
    myfile.close()
    for i in range(7):
        s[i] = s[i].split(',')

    s = s[:-1]

    for i in range(7):
        for j in range(len(s[i])):
            s[i][j] = float(s[i][j])

    v_mag_earth = [[0,0,0] for i in range(7)]
    v_grav_earth = [[0,0,0] for i in range(7)]
    v_mag_mobile = [[0,0,0] for i in range(7)]
    v_grav_mobile = [[0,0,0] for i in range(7)]
    for i in range(7):
        v_mag_earth[i] = array([s[i][0], s[i][1], s[i][2]])
        v_grav_earth[i] = array([s[i][3], s[i][4], s[i][5]])
        v_mag_mobile[i] = array([s[i][6], s[i][7], s[i][8]])
        v_grav_mobile[i] = array([s[i][9], s[i][10], s[i][11]])


    #for i in range(7):
    #    QuEST01(v_mag_earth[i],v_grav_earth[i],v_mag_mobile[i],v_grav_mobile[i],0.9)
    Quar = [0,0,0,0]
    Quar = QuEST01([0,-0.913,0.411],[0,0,-1],[0,-0.911,0.411],[0,0,-1],0.5)
    theta = 0.001
    e = [0.01,0.01,0.01]
    theta = 2.0*acos(Quar[3])
    e[0] = Quar[0]/sin(theta/2)
    e[1] = Quar[1]/sin(theta/2)
    e[2] = Quar[2]/sin(theta/2)
    print e,theta*180/pi
