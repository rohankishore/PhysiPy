from .constants import *


class Gravitation:
    def __init__(self):
        super().__init__()

    def gr(m1, m2, r):
        rsq = r * r
        result = G * ((m1 * m2) / rsq)
        return result

    def G_Potential(M, r):
        v = (-G * M) / r
        return v

    def g_in_depth(depth):
        return g * (1 - (depth / 6400))

    def axial_velocity(area_swept, time):
        return area_swept / time

    def gravitational_force(mass1, mass2, distance):
        return G * mass1 * mass2 / distance ** 2

    @staticmethod
    def gravitational_potential_energy(mass1, mass2, distance):
        return -(G * mass1 * mass2) / distance

    @staticmethod
    def gravitational_field_strength(mass, distance):
        return G * mass / distance ** 2

    @staticmethod
    def escape_velocity(mass, radius):
        return math.sqrt(2 * G * mass / radius)

    @staticmethod
    def orbital_velocity(mass, radius):
        return math.sqrt(G * mass / radius)

    @staticmethod
    def period_of_orbit(mass, radius):
        return 2 * math.pi * math.sqrt(radius ** 3 / (G * mass))

    @staticmethod
    def gravitational_potential(mass, distance):
        return -(G * mass) / distance

    @staticmethod
    def weight(mass, acceleration_due_to_gravity):
        return mass * acceleration_due_to_gravity

    @staticmethod
    def gravitational_acceleration(mass, distance):
        return G * mass / distance ** 2

    @staticmethod
    def keplers_third_law(period, semi_major_axis):
        return (period ** 2) / (semi_major_axis ** 3)
