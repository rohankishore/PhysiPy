from .constants import * 
import math


class Thermodynamics:
    @staticmethod
    def temperature_conversion_celsius_to_kelvin(celsius):
        return celsius + 273.15

    @staticmethod
    def temperature_conversion_kelvin_to_celsius(kelvin):
        return kelvin - 273.15

    @staticmethod
    def ideal_gas_law(pressure, volume, temperature):
        return (pressure * volume) / (boltzmann_constant * temperature)

    @staticmethod
    def thermal_expansion_coefficient(change_in_length, original_length, change_in_temperature):
        return change_in_length / (original_length * change_in_temperature)

    @staticmethod
    def heat_transfer_conduction(thermal_conductivity, area, temperature_difference, thickness):
        return (thermal_conductivity * area * temperature_difference) / thickness

    @staticmethod
    def heat_transfer_convection(heat_transfer_coefficient, area, temperature_difference):
        return heat_transfer_coefficient * area * temperature_difference

    @staticmethod
    def heat_transfer_radiation(emissivity, stefan_boltzmann_constant, area, temperature1, temperature2):
        return emissivity * stefan_boltzmann_constant * area * (temperature1 ** 4 - temperature2 ** 4)

    @staticmethod
    def first_law_thermodynamics_bool(heat_added, work_done, change_in_internal_energy):
        return heat_added + work_done == change_in_internal_energy

    @staticmethod
    def first_law_thermodynamics(heat_added, work_done):
        return heat_added + work_done

    @staticmethod
    def efficiency_carnot(temperature_hot, temperature_cold):
        return 1 - temperature_cold / temperature_hot

    @staticmethod
    def efficiency_heat_engine(heat_added, heat_rejected):
        return (heat_added - heat_rejected) / heat_added

    @staticmethod
    def entropy_change(heat_transfer, temperature):
        return heat_transfer / temperature

    @staticmethod
    def entropy_change_irreversible(heat_transfer, temperature_reservoir):
        return heat_transfer / temperature_reservoir

    @staticmethod
    def entropy_change_adiabatic(reversible_entropy_change, heat_transfer, temperature_reservoir):
        return reversible_entropy_change - heat_transfer / temperature_reservoir

    @staticmethod
    def entropy_change_phase(heat_transfer, temperature):
        return heat_transfer / temperature

    @staticmethod
    def work_done_by_ideal_gas(pressure, volume_initial, volume_final):
        return pressure * (volume_final - volume_initial)

    @staticmethod
    def root_mean_square_speed(molar_mass, temperature):
        return math.sqrt((3 * boltzmann_constant * temperature) / molar_mass)

    @staticmethod
    def average_kinetic_energy(temperature):
        return (3 / 2) * boltzmann_constant * temperature

    @staticmethod
    def average_kinetic_energy_with_molar_mass(molar_mass, temperature):
        return (3 / 2) * boltzmann_constant * molar_mass * temperature

    @staticmethod
    def speed_of_sound(molar_mass, temperature, gamma):
        return math.sqrt(gamma * boltzmann_constant * temperature / molar_mass)

    @staticmethod
    def specific_heat_capacity(mass, temperature_change, heat_transfer):
        return heat_transfer / (mass * temperature_change)

    @staticmethod
    def latent_heat(mass, heat_transfer):
        return heat_transfer / mass
