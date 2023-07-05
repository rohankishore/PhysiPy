import numpy as np
import math

g = 9.8
G = 6.67408 * (math.pow(10, -11))
speed_of_light = 3 * 10 ^ 8

def weight(mass):
    result = mass * 9.8
    return result

def density(mass, volume):
    result = mass/volume
    return result

def volume(L, B, H):
    print("Volume:", L * B * H)

def area(L, B):
    print("Area:", L * B)

def error_muldiv(a, b, c, d):
    x = ((a / c) + (b / d))
    return x

def error_addsub(a, b):
    p = (a + b)
    return p

def percentage_error(M, E):
    print("Percentage Error:", (M / E) * 100, "%")

def absolute_error(a, n):
    if len(a) == n:
        A = np.array(a)
        S = np.sum(A)
        E = S / n
        return E
    else:
        return ValueError

def meanabsolute_error(a, n):
    if len(a) == n:
        b = []
        for i in a:
            b.append(abs(i))
        A = np.array(b)
        S = np.sum(A)
        M = S / n
        return M
    else:
        print("error try again")

def speed(distance, time):
    result = distance/time
    return result

def velocity(displacement, time):
    return displacement/time

def displacement(velocity, time):
    return velocity*time

def force(mass, acceleration):
    return mass * acceleration

def momentum(mass, velocity):
    return mass * velocity

def distance(speed, time):
    return speed*time

def acceleration(u,v, dt):
    return (v - u)/dt

def time(distance, speed):
    return distance / speed


def G(m1, m2, r):
    rsq = r*r
    result = G * ((m1*m2)/rsq)
    return result

def G_Potential(M, r):
    v = (-G * M) / r
    return v

def g_in_depth(depth):
    res = g * (1-(depth/6400))
    return res

def axial_velocity(area_swept, time):
    return area_swept/time
