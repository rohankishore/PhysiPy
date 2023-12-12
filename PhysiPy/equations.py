from .constants import g


def weight(mass):
    return mass * g


def density(m, v):
    return m / v


def volume(l, w, h):
    return l * w * h


def area(l, w):
    return l * w
