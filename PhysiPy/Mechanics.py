from .constants import pi
import math


class Mechanics:

    def __init__(self, initial_velocity=1, acceleration=1, time=1, final_velocity=1,
                 mass=1, velocity=1, force=1, displacement=1, angle=0, work=1,
                 radius=1, period=1, lever_arm=1, angular_displacement=1,
                 angular_velocity=1):
        self.initial_velocity = initial_velocity
        self.acceleration = acceleration
        self.time = time
        self.final_velocity = final_velocity
        self.mass = mass
        self.velocity = velocity
        self.force = force
        self.displacement = displacement
        self.angle = angle
        self.work = work
        self.radius = radius
        self.period = period
        self.lever_arm = lever_arm
        self.angular_displacement = angular_displacement
        self.angular_velocity = angular_velocity

    def velocity(self):
        return self.initial_velocity + self.acceleration * self.time

    def displacement(self):
        return self.initial_velocity * self.time + 0.5 * self.acceleration * \
               self.time ** 2

    def acceleration(self):
        return (self.final_velocity - self.initial_velocity) / self.time

    def uniform_accelerated_motion(self):
        displacement = self.initial_velocity * self.time + 0.5 * self.acceleration * \
                       self.time ** 2
        final_velocity = self.initial_velocity + self.acceleration * self.time
        return displacement, final_velocity

    def force(self):
        return self.mass * self.acceleration

    def work(self):
        return self.force * self.displacement * math.cos(math.radians(self.angle))

    def kinetic_energy(self):
        return 0.5 * self.mass * self.velocity ** 2

    def potential_energy(self, height, gravitational_field_strength):
        return self.mass * height * gravitational_field_strength

    def power(self):
        return self.work / self.time

    def momentum(self):
        return self.mass * self.velocity

    def impulse(self):
        return self.force * self.time

    def circular_velocity(self):
        return 2 * pi * self.radius / self.period

    def centripetal_acceleration(self):
        return self.velocity ** 2 / self.radius

    def torque(self):
        return self.force * self.lever_arm

    def angular_velocity(self):
        return self.angular_displacement / self.time

    def angular_acceleration(self):
        return self.angular_velocity / self.time
