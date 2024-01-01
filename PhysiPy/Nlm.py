# from .constants import *


class Nlm:
    def __init__(self, mass=1, acceleration=1, velocity=1, massofbullet=1,
                 massofgun=1, initialvelocity=1):
        self.mass = mass
        self.acceleration = acceleration
        self.velocity = velocity
        self.massofbullet = massofbullet
        self.massofgun = massofgun
        self.initialvelocity = initialvelocity

    def force(self):
        return self.mass * self.acceleration

    def momentum(self):
        return self.mass * self.velocity

    def recoil_velocity(self):
        return (self.massofbullet * self.initialvelocity) / self.massofgun
