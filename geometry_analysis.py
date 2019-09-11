from coordinates import *
from correlation import *
import math
import numpy as np
import quaternion as quat

def cosrule(a, b, c, deg = False):
    alpha = math.acos((b**2 + c**2 - a**2) / (2*b*c))
    if deg:
        alpha *= (180 / math.pi)
    return alpha

def angle3p(a, b, c, deg = False):
    s_ab = np.linalg.norm(a - b)
    s_bc = np.linalg.norm(b - c)
    s_ca = np.linalg.norm(c - a)

    alpha = cosrule(s_bc, s_ab, s_ca, deg)
    return alpha

def rotate(v, axis, theta):
    theta *= (math.pi / 180)
    vector = np.array([0.] + v)
    rot_axis = np.array([0.] + axis)
    axis_angle = (theta*0.5) * rot_axis/np.linalg.norm(rot_axis)

    vec = quat.quaternion(*v)
    qlog = quat.quaternion(*axis_angle)
    q = np.exp(qlog)
    v_prime = q * vec * np.conjugate(q)

    return v_prime.imag
