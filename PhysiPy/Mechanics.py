from .constants import *
import math


class Mechanics:
    def __init__(self):
        super().__init__()

    @staticmethod
    def velocity(initial_velocity, acceleration, time):
        return initial_velocity + acceleration * time

    @staticmethod
    def displacement(initial_velocity, acceleration, time):
        return initial_velocity * time + 0.5 * acceleration * time ** 2

    @staticmethod
    def acceleration(final_velocity, initial_velocity, time):
        return (final_velocity - initial_velocity) / time

    @staticmethod
    def uniform_accelerated_motion(initial_velocity, time, acceleration):
        displacement = initial_velocity * time + 0.5 * acceleration * time ** 2
        final_velocity = initial_velocity + acceleration * time
        return displacement, final_velocity

    @staticmethod
    def force(mass, acceleration):
        return mass * acceleration

    @staticmethod
    def work(force, displacement, angle=0):
        return force * displacement * math.cos(math.radians(angle))

    @staticmethod
    def kinetic_energy(mass, velocity):
        return 0.5 * mass * velocity ** 2

    @staticmethod
    def potential_energy(mass, height, gravitational_field_strength):
        return mass * height * gravitational_field_strength

    @staticmethod
    def power(work, time):
        return work / time

    @staticmethod
    def momentum(mass, velocity):
        return mass * velocity

    @staticmethod
    def impulse(force, time):
        return force * time

    @staticmethod
    def circular_velocity(radius, period):
        return 2 * math.pi * radius / period

    @staticmethod
    def centripetal_acceleration(radius, velocity):
        return velocity ** 2 / radius

    @staticmethod
    def torque(force, lever_arm):
        return force * lever_arm

    @staticmethod
    def angular_velocity(angular_displacement, time):
        return angular_displacement / time

    @staticmethod
    def angular_acceleration(angular_velocity, time):
        return angular_velocity / time
