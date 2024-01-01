from .constants import boltzmann_constant, stefan_boltzmann_constant
import math


class Thermodynamics:

    def __init__(self, celsius=1, kelvin=1, heat_added=1, work_done=1,
                 change_in_internal_energy=1, pressure=1, volume=1, temperature=1,
                 change_in_length=1, original_length=1, change_in_temperature=1,
                 thermal_conductivity=1, area=1, temperature_difference=1, thickness=1,
                 molar_mass=1, gamma=1, mass=1, heat_transfer=1, temperature_change=1,
                 temperature_hot=1, temperature_cold=1, heat_rejected=1,
                 volume_initial=1, volume_final=1, heat_transfer_coefficient=1,
                 emissivity=1, temperature1=1, temperature2=1, temperature_reservoir=1,
                 reversible_entropy_change=1):
        self.celsius = celsius
        self.kelvin = kelvin
        self.heat_added = heat_added
        self.heat_rejected = heat_rejected
        self.work_done = work_done
        self.change_in_internal_energy = change_in_internal_energy
        self.pressure = pressure
        self.volume = volume
        self.temperature = temperature
        self.change_in_length = change_in_length
        self.original_length = original_length
        self.change_in_temperature = change_in_temperature
        self.thermal_conductivity = thermal_conductivity
        self.area = area
        self.temperature_difference = temperature_difference
        self.thickness = thickness
        self.molar_mass = molar_mass
        self.gamma = gamma
        self.mass = mass
        self.heat_transfer = heat_transfer
        self.temperature_change = temperature_change
        self.temperature_hot = temperature_hot
        self.temperature_cold = temperature_cold
        self.volume_initial = volume_initial
        self.volume_final = volume_final
        self.heat_transfer_coefficient = heat_transfer_coefficient
        self.emissivity = emissivity
        self.temperature1 = temperature1
        self.temperature2 = temperature2
        self.temperature_reservoir = temperature_reservoir
        self.reversible_entropy_change = reversible_entropy_change

    def convert_celsius_to_kelvin(self) -> float:
        """Convert temperature from Celsius to Kelvin."""
        return self.celsius + 273.15

    def convert_kelvin_to_celsius(self) -> float:
        """Convert temperature from Kelvin to Celsius."""
        return self.kelvin - 273.15

    def ideal_gas_law(self) -> float:
        return (self.pressure * self.volume) / (boltzmann_constant * self.temperature)

    def thermal_expansion_coefficient(self) -> float:
        return self.change_in_length / (self.original_length *
                                        self.change_in_temperature)

    def heat_transfer_conduction(self) -> float:
        return (self.thermal_conductivity * self.area * self.temperature_difference) / \
               self.thickness

    def heat_transfer_convection(self) -> float:
        return self.heat_transfer_coefficient * self.area * self.temperature_difference

    def heat_transfer_radiation(self) -> float:
        return self.emissivity * stefan_boltzmann_constant * self.area * \
               (self.temperature1 ** 4 - self.temperature2 ** 4)

    def is_first_law_satisfied(self) -> bool:
        return self.first_law_thermodynamics() == self.change_in_internal_energy

    def first_law_thermodynamics(self) -> float:
        return self.heat_added + self.work_done

    def efficiency_carnot(self) -> float:
        return 1 - self.temperature_cold / self.temperature_hot

    def efficiency_heat_engine(self) -> float:
        return (self.heat_added - self.heat_rejected) / self.heat_added

    def entropy_change(self) -> float:
        return self.heat_transfer / self.temperature

    def entropy_change_irreversible(self) -> float:
        return self.heat_transfer / self.temperature_reservoir

    def entropy_change_adiabatic(self) -> float:
        return self.reversible_entropy_change - self.heat_transfer / \
               self.temperature_reservoir

    def entropy_change_phase(self) -> float:
        return self.heat_transfer / self.temperature

    def work_done_by_ideal_gas(self) -> float:
        return self.pressure * (self.volume_final - self.volume_initial)

    def root_mean_square_speed(self) -> float:
        return math.sqrt((3 * boltzmann_constant * self.temperature) / self.molar_mass)

    def average_kinetic_energy(self) -> float:
        return (3 / 2) * boltzmann_constant * self.temperature

    def average_kinetic_energy_with_molar_mass(self) -> float:
        return (3 / 2) * boltzmann_constant * self.molar_mass * self.temperature

    def speed_of_sound(self) -> float:
        return math.sqrt(self.gamma * boltzmann_constant * self.temperature /
                         self.molar_mass)

    def specific_heat_capacity(self) -> float:
        return self.heat_transfer / (self.mass * self.temperature_change)

    def latent_heat(self) -> float:
        return self.heat_transfer / self.mass
