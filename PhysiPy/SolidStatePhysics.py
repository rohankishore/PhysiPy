from .constants import *
import math
from .Mass import Mass
from .Charge import Charge


class SolidStatePhysics:
    @staticmethod
    def electrical_resistivity(resistance, cross_sectional_area, length):
        return resistance * (cross_sectional_area / length)

    @staticmethod
    def hall_effect(hall_voltage, current, magnetic_field, thickness):
        return (hall_voltage * thickness) / (current * magnetic_field)

    @staticmethod
    def electron_mobility(hall_coefficient, electron_charge, electron_density):
        return hall_coefficient / (electron_charge * electron_density)

    @staticmethod
    def fermi_energy(electron_density, work_function):
        return (electron_density * (plancks_constant ** 2)) / (2 * Mass.electron) + work_function

    @staticmethod
    def band_gap(fermi_energy, temperature):
        return 2 * boltzmann_constant * temperature * math.log(2) - fermi_energy

    @staticmethod
    def energy_band(kinetic_energy, potential_energy):
        return kinetic_energy + potential_energy

    @staticmethod
    def energy_band_conduction(fermi_energy, intrinsic_fermi_level):
        return fermi_energy - intrinsic_fermi_level

    @staticmethod
    def intrinsic_carrier_concentration(donor_density, acceptor_density, fermi_energy, intrinsic_fermi_level,
                                        temperature):
        return ((donor_density * acceptor_density) ** 0.5) * math.exp(
            -(fermi_energy - intrinsic_fermi_level) / (2 * boltzmann_constant * temperature))

    @staticmethod
    def depletion_layer_width(fermi_energy, potential_difference, electron_density):
        return math.sqrt((2 * permittivity_of_free_space * potential_difference) / (
                Charge.electron * electron_density)) / (2 * fermi_energy)

    @staticmethod
    def solar_cell_efficiency(power_output, light_power_input):
        return (power_output / light_power_input) * 100

    @staticmethod
    def energy_density(material_density, specific_heat_capacity, temperature_change):
        return material_density * specific_heat_capacity * temperature_change

    @staticmethod
    def youngs_modulus(stress, strain):
        return stress / strain

    @staticmethod
    def poisson_ratio(transverse_strain, longitudinal_strain):
        return -transverse_strain / longitudinal_strain

    @staticmethod
    def thermal_conductivity(heat_flux, temperature_gradient, thickness, cross_sectional_area):
        return heat_flux * thickness / (temperature_gradient * cross_sectional_area)

    @staticmethod
    def superconductivity_critical_temperature(critical_field):
        return (critical_field ** 2) * (
                    permittivity_of_free_space / (2 * Mass.electron * Charge.electron))

    @staticmethod
    def superconductivity_coherence_length(hamiltonian_operator, fermi_energy):
        return (plancks_constant / (2 * math.pi)) * (
                (3 * permittivity_of_free_space * (hamiltonian_operator ** 2)) / (
                    2 * Mass.electron * fermi_energy)) ** 0.5
