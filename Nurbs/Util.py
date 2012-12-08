import math

dependencies = '''This module requires:
	Numeric Python
'''

try:
    import numpy as np
except ImportError, value:
    print dependencies
    raise

def uniformknots(cntrlpts, degree):
    knots = np.zeros(degree + cntrlpts + 1, np.float)
    knots[cntrlpts:] = 1.
    knots[degree+1:cntrlpts] = np.arange(1., cntrlpts-degree)* (1./(cntrlpts - degree))
    return knots

def translate(txyz):
    ret = np.identity(4).astype(np.float)
    ret[0:len(txyz), 3] = txyz
    return ret

def scale(sxyz):
    ret = np.identity(4).astype(np.float)
    s = np.ones(3, np.float)
    s[0:len(sxyz)] = sxyz
    ret[0,0] = s[0]
    ret[1,1] = s[1]
    ret[2,2] = s[2]
    return ret

def deg2rad(angle):
    return math.pi * angle / 180.

def rad2deg(angle):
    return angle * 180./math.pi

def rotx(angle):
    ret = np.identity(4).astype(np.float)
    ret[1,1] = ret[2,2] = np.cos(angle)
    sn = np.sin(angle)
    ret[1,2] = -sn
    ret[2,1] = sn
    return ret

def roty(angle):
    ret = np.identity(4).astype(np.float)
    ret[0,0] = ret[2,2] = np.cos(angle)
    sn = np.sin(angle)
    ret[0,2] = sn
    ret[2,0] = -sn
    return ret

def rotz(angle):
    ret = np.identity(4).astype(np.float)
    ret[0,0] = ret[1,1] = np.cos(angle)
    sn = np.sin(angle)
    ret[0,1] = -sn
    ret[1,0] = sn
    return ret

def rotxyz(angles):
    ret = np.identity(4).astype(np.float)
    
    ret[1,1] = ret[2,2] = np.cos(angles[0])
    sn = np.sin(angles[0])
    ret[1,2] = -sn
    ret[2,1] = sn
    
    cs = np.cos(angles[1])
    ret[0,0] = cs
    ret[2,2] += cs
    sn = np.sin(angles[1])
    ret[0,2] = sn
    ret[2,0] = -sn

    cs = np.cos(angles[2])
    ret[0,0] += cs
    ret[1,1] += cs
    sn = np.sin(angles[2])
    ret[0,1] = -sn
    ret[1,0] = sn
    return ret

