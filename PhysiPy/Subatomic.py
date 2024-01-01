# from .constants import *
import math


class Subatomic:
    def __init__(self, mass=1, mass_parent=1, mass_daughters=1, momentum=1,
                 atomic_number=1, n=1, initial_amount=1, decay_constant=1, time=1):
        self.mass = mass
        self.mass_parent = mass_parent
        self.mass_daughters = mass_daughters
        self.momentum = momentum
        self.atomic_number = atomic_number
        self.n = n
        self.initial_amount = initial_amount
        self.decay_constant = decay_constant
        self.time = time

    def mass_energy_equivalence(self):
        return self.mass * 3e8 ** 2

    def binding_energy(self):
        # Energy in electron volts (eV)
        return (self.mass_parent - sum(self.mass_daughters)) * 9e16

    def de_broglie_wavelength(self) -> float:
        return 6.63e-34 / (self.momentum * self.mass)

    def bohr_radius(self):
        # Distance in angstroms (Å)
        return 0.529 / self.atomic_number

    def energy_level_hydrogen(self):
        # Energy in electron volts (eV)
        return -13.6 / (self.n ** 2)

    def radioactive_decay(self):
        return self.initial_amount * math.exp(-self.decay_constant * self.time)

    def half_life(self):
        return math.log(2) / self.decay_constant
