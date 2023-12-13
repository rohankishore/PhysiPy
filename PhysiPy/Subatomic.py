from .constants import *
import math


class Subatomic:
    def __init__(self):
        super().__init__()

    @staticmethod
    def mass_energy_equivalence(mass):
        return mass * 3e8 ** 2

    @staticmethod
    def binding_energy(mass_parent, mass_daughters):
        # Energy in electron volts (eV)
        return (mass_parent - sum(mass_daughters)) * 9e16

    @staticmethod
    def de_broglie_wavelength(momentum, mass):
        return 6.63e-34 / (momentum * mass)

    @staticmethod
    def bohr_radius(atomic_number):
        # Distance in angstroms (Å)
        return 0.529 / atomic_number

    @staticmethod
    def energy_level_hydrogen(n):
        # Energy in electron volts (eV)
        return -13.6 / (n ** 2)

    @staticmethod
    def radioactive_decay(initial_amount, decay_constant, time):
        return initial_amount * math.exp(-decay_constant * time)

    @staticmethod
    def half_life(decay_constant):
        return math.log(2) / decay_constant
