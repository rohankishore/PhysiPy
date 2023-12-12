from .constants import *
import math


class Subatomic:
    def __init__(self):
        super().__init__()

    @staticmethod
    def mass_energy_equivalence(mass):
        return mass * (3e8) ** 2

    @staticmethod
    def binding_energy(mass_parent, mass_daughters):
        return (mass_parent - sum(mass_daughters)) * (9e16)  # Energy in electron volts (eV)

    @staticmethod
    def de_broglie_wavelength(momentum, mass):
        return (6.63e-34) / (momentum * mass)

    @staticmethod
    def bohr_radius(atomic_number):
        return 0.529 / atomic_number  # Distance in angstroms (Å)

    @staticmethod
    def energy_level_hydrogen(n):
        return (-13.6) / (n ** 2)  # Energy in electron volts (eV)

    @staticmethod
    def radioactive_decay(initial_amount, decay_constant, time):
        return initial_amount * math.exp(-decay_constant * time)

    @staticmethod
    def half_life(decay_constant):
        return math.log(2) / decay_constant

    class Mass:
        electron = 9.1 * (math.pow(10, -31))
        proton = 1.67 * (math.pow(10, -27))
        neutron = 1.67 * (math.pow(10, -27))

    class Charge:
        electron = 1.602 * (math.pow(10, -19))
