from .constants import plancks_constant, Coulombs_constant, boltzmann_constant
import math
from .Mass import Mass
from .Charge import Charge


class SolidStatePhysics:

    def __init__(self, resistance=1, cross_sectional_area=1, length=1, hall_voltage=1,
                 current=1, magnetic_field=1, thickness=1, hall_coefficient=1,
                 electron_charge=1, electron_density=1, work_function=1, fermi_energy=1,
                 temperature=1, kinetic_energy=1, potential_energy=1,
                 intrinsic_fermi_level=1, donor_density=1, acceptor_density=1,
                 potential_difference=1, power_output=1, light_power_input=1,
                 material_density=1, specific_heat_capacity=1, temperature_change=1,
                 stress=1, strain=1, transverse_strain=1, longitudinal_strain=1,
                 heat_flux=1, temperature_gradient=1, critical_field=1,
                 hamiltonian_operator=1):
        self.r = resistance
        self.cross_sectional_area = cross_sectional_area
        self.length = length
        self.hall_voltage = hall_voltage
        self.C = current
        self.magnetic_field = magnetic_field
        self.thickness = thickness
        self.hall_coefficient = hall_coefficient
        self.electron_charge = electron_charge
        self.electron_density = electron_density
        self.work_function = work_function
        self.fermi_energy = fermi_energy
        self.temperature = temperature
        self.kinetic_energy = kinetic_energy
        self.potential_energy = potential_energy
        self.intrinsic_fermi_level = intrinsic_fermi_level
        self.donor_density = donor_density
        self.acceptor_density = acceptor_density
        self.potential_difference = potential_difference
        self.power_output = power_output
        self.light_power_input = light_power_input
        self.material_density = material_density
        self.specific_heat_capacity = specific_heat_capacity
        self.temperature_change = temperature_change
        self.stress = stress
        self.strain = strain
        self.transverse_strain = transverse_strain
        self.longitudinal_strain = longitudinal_strain
        self.heat_flux = heat_flux
        self.temperature_gradient = temperature_gradient
        self.critical_field = critical_field
        self.hamiltonian_operator = hamiltonian_operator

    def electrical_resistivity(self):
        return self.r * (self.cross_sectional_area / self.length)

    def hall_effect(self):
        return (self.hall_voltage * self.thickness) / (self.C * self.magnetic_field)

    def electron_mobility(self):
        return self.hall_coefficient / (self.electron_charge * self.electron_density)

    def fermi_energy(self):
        return (self.electron_density * (plancks_constant ** 2)) / \
               (2 * Mass.electron) + self.work_function

    def band_gap(self):
        return 2 * boltzmann_constant * self.temperature * math.log(2) - \
               self.fermi_energy

    def energy_band(self):
        return self.kinetic_energy + self.potential_energy

    def energy_band_conduction(self):
        return self.fermi_energy - self.intrinsic_fermi_level

    def intrinsic_carrier_concentration(self):
        return ((self.donor_density * self.acceptor_density) ** 0.5) * \
               math.exp(-(self.fermi_energy - self.intrinsic_fermi_level) /
                        (2 * boltzmann_constant * self.temperature))

    def depletion_layer_width(self):
        return math.sqrt((2 * Coulombs_constant * self.potential_difference) / (
                Charge.electron * self.electron_density)) / (2 * self.fermi_energy)

    def solar_cell_efficiency(self):
        return (self.power_output / self.light_power_input) * 100

    def energy_density(self):
        return self.material_density * self.specific_heat_capacity * \
               self.temperature_change

    def youngs_modulus(self):
        return self.stress / self.strain

    def poisson_ratio(self):
        return -self.transverse_strain / self.longitudinal_strain

    def thermal_conductivity(self):
        return self.heat_flux * self.thickness / (self.temperature_gradient *
                                                  self.cross_sectional_area)

    def superconductivity_critical_temperature(self):
        return (self.critical_field ** 2) * (
                    Coulombs_constant / (2 * Mass.electron * Charge.electron))

    def superconductivity_coherence_length(self):
        return (plancks_constant / (2 * math.pi)) * ((3 * Coulombs_constant *
                (self.hamiltonian_operator ** 2)) / (2 * Mass.electron *
                                                     self.fermi_energy)) ** 0.5
