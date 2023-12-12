from .constants import *

class Nlm:
    def __init__(self):
        super().__init__()

    @staticmethod
    def force(mass, acceleration):
        return mass * acceleration

    @staticmethod
    def momentum(mass, velocity):
        return mass * velocity

    @staticmethod
    def recoil_velocity(massofbullet, massofgun, initalvelocity):
        return (massofbullet * initalvelocity) / massofgun
