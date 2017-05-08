from numpy import *
from numpy.linalg import norm, inv
from Euler_from_Quaternion01 import Quat

def give_rotated(Quaternion, initial_v):
    initial_v = array(initial_v)
    Quaternion = Quaternion/norm(Quaternion)
    Q = Quat(Quaternion)
    rmat = Q._quat2transform()
    rotated = rmat.dot(initial_v.T)
    print rotated #,"\n", rmat
give_rotated([0,0,0.259,0.966], [41.5, -1.3, -16])
