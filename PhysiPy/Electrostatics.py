from .constants import *


class Electrostatics:
    @staticmethod
    def electric_force(charge1, charge2, distance):
        return Coulombs_constant * charge1 * charge2 / distance ** 2

    @staticmethod
    def electric_field(charge, distance):
        return Coulombs_constant * charge / distance ** 2

    @staticmethod
    def electric_potential(charge, distance):
        return Coulombs_constant * charge / distance

    @staticmethod
    def capacitance(charge, voltage):
        return charge / voltage

    @staticmethod
    def electric_current(charge, time):
        return charge / time

    @staticmethod
    def resistance(voltage, current):
        return voltage / current

    @staticmethod
    def ohms_law(voltage, current):
        return voltage / current

    @staticmethod
    def coulombs_law(electric_field, charge):
        return electric_field * charge

    @staticmethod
    def gauss_law(charge_enclosed, epsilon_0):
        return charge_enclosed / epsilon_0

    @staticmethod
    def faradays_law(emf, time):
        return emf * time

    @staticmethod
    def magnetic_field(force, charge, velocity):
        return force / (charge * velocity)

    @staticmethod
    def lorentz_force(charge, velocity, magnetic_field):
        return charge * velocity * magnetic_field

    @staticmethod
    def hall_voltage(magnetic_field, current, charge_carrier_density, thickness):
        return (magnetic_field * current) / (charge_carrier_density * thickness)

    @staticmethod
    def drift_velocity(current, charge, charge_carrier_density, cross_sectional_area):
        return current / (charge * charge_carrier_density * cross_sectional_area)

    @staticmethod
    def resistivity(resistance, cross_sectional_area, length):
        return resistance * (cross_sectional_area / length)

    def series_resistance(*resistances):
        return sum(resistances)

    @staticmethod
    def parallel_resistance(*resistances):
        return 1 / sum(1 / r for r in resistances)
