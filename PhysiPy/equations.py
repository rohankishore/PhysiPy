from .constants import g


def weight(mass):
    return mass * g


def density(mass, vol):
    return mass / vol


def volume(length, width, height):
    return length * width * height


def area(length, width):
    return length * width
