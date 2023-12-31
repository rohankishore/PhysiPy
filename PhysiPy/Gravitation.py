from .constants import g, G, pi
import math


class Gravitation:

    def __init__(self, m1=1, m2=1, r=1, M=1, depth=1, area_swept=1, time=1, mass1=1,
                 mass2=1, mass=1, distance=1, radius=1,
                 acceleration_due_to_gravity=1, period=1, semi_major_axis=1):
        self.m1 = m1
        self.m2 = m2
        self.r = r
        self.M = M
        self.depth = depth
        self.area_swept = area_swept
        self.time = time
        self.mass1 = mass1
        self.mass2 = mass2
        self.mass = mass
        self.distance = distance
        self.radius = radius
        self.acceleration_due_to_gravity = acceleration_due_to_gravity
        self.period = period
        self.semi_major_axis = semi_major_axis

    def gravity(self):
        rsq = self.r * self.r
        return G * ((self.m1 * self.m2) / rsq)

    def G_Potential(self):
        return (-G * self.M) / self.r

    def g_in_depth(self):
        return g * (1 - (self.depth / 6400))

    def axial_velocity(self):
        return self.area_swept / self.time

    def gravitational_force(self):
        return G * self.mass1 * self.mass2 / self.distance ** 2

    def gravitational_potential_energy(self):
        return -(G * self.mass1 * self.mass2) / self.distance

    def gravitational_field_strength(self):
        return G * self.mass / self.distance ** 2

    def escape_velocity(self):
        return math.sqrt(2 * G * self.mass / self.radius)

    def orbital_velocity(self):
        return math.sqrt(G * self.mass / self.radius)

    def period_of_orbit(self):
        return 2 * pi * math.sqrt(self.radius ** 3 / (G * self.mass))

    def gravitational_potential(self):
        return -(G * self.mass) / self.distance

    def weight(self):
        return self.mass * self.acceleration_due_to_gravity

    def gravitational_acceleration(self):
        return G * self.mass / self.distance ** 2

    def keplers_third_law(self):
        return (self.period ** 2) / (self.semi_major_axis ** 3)
