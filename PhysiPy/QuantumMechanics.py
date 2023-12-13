from .constants import *
import math
from .Mass import Mass
from .Charge import Charge


class QuantumMechanics:
    @staticmethod
    def uncertainty_principle(position_uncertainty, momentum_uncertainty):
        return position_uncertainty * momentum_uncertainty >= plancks_constant / (4 * math.pi)

    @staticmethod
    def schrodingers_equation(wave_function, hamiltonian_operator, energy):
        return hamiltonian_operator * wave_function == energy * wave_function

    @staticmethod
    def probability_density(wave_function):
        return wave_function.conjugate() * wave_function

    @staticmethod
    def tunneling_probability(energy_barrier, particle_mass, particle_energy):
        return math.exp((-2 * math.sqrt(2 * particle_mass * energy_barrier) / plancks_constant) *
                        math.sqrt(particle_energy - energy_barrier))

    @staticmethod
    def black_body_radiation_intensity(frequency, temperature):
        return (2 * plancks_constant * frequency ** 3) / (
                speed_of_light ** 2 * (
                math.exp((plancks_constant * frequency) / (boltzmann_constant * temperature)) - 1))

    @staticmethod
    def black_body_radiation_power(temperature):
        return (2 * math.pi ** 5 * boltzmann_constant ** 4 * temperature ** 4) / (15 * plancks_constant ** 3 *
                                                                                  speed_of_light ** 2)

    @staticmethod
    def fermi_dirac_distribution(energy, fermi_level, temperature):
        return 1 / (math.exp((energy - fermi_level) / (boltzmann_constant * temperature)) + 1)

    @staticmethod
    def bose_einstein_distribution(energy, chemical_potential, temperature):
        return 1 / (math.exp((energy - chemical_potential) / (boltzmann_constant * temperature)) - 1)

    @staticmethod
    def wave_particle_duality(wavelength, momentum):
        return wavelength * momentum == plancks_constant

    @staticmethod
    def heisenberg_uncertainty_energy_lifetime(energy, lifetime):
        return energy * lifetime >= plancks_constant / (4 * math.pi)

    @staticmethod
    def heisenberg_uncertainty_position_momentum(position_uncertainty, momentum_uncertainty):
        return position_uncertainty * momentum_uncertainty >= plancks_constant / 2

    @staticmethod
    def atomic_orbital_radius(principal_quantum_number):
        return (4 * math.pi * permittivity_of_free_space * (plancks_constant ** 2) * (
                principal_quantum_number ** 2)) / (Mass.electron * (elementary_charge ** 2))

    @staticmethod
    def fine_structure_constant():
        return (elementary_charge ** 2) / (4 * math.pi * permittivity_of_free_space * plancks_constant * speed_of_light)

    @staticmethod
    def de_broglie_wavelength_photon(frequency):
        return speed_of_light / frequency

    @staticmethod
    def de_broglie_wavelength_particle(momentum):
        return plancks_constant / momentum

    @staticmethod
    def plank_distribution_radiation(energy, temperature):
        return (2 * plancks_constant * (energy ** 3)) / (
                speed_of_light ** 2 * (math.exp((plancks_constant * energy) / (boltzmann_constant * temperature)) - 1))

    @staticmethod
    def plank_distribution(energy, temperature):
        return (2 * energy ** 2) / (
                (plancks_constant ** 3) * (speed_of_light ** 2) *
                (math.exp(energy / (boltzmann_constant * temperature)) - 1))

    @staticmethod
    def plank_law_intensity(frequency, temperature):
        return (2 * plancks_constant * frequency ** 3) / (
                speed_of_light ** 2 * (
                math.exp((plancks_constant * frequency) / (boltzmann_constant * temperature)) - 1))

    @staticmethod
    def plank_law_power(temperature):
        return (2 * math.pi * (boltzmann_constant ** 4) * temperature ** 4) / (15 * (plancks_constant ** 3) *
                                                                               (speed_of_light ** 2))

    @staticmethod
    def einstein_light_quanta(energy, frequency):
        return energy == plancks_constant * frequency

    @staticmethod
    def einstein_light_intensity(number_of_quanta, time):
        return number_of_quanta / time

    @staticmethod
    def atomic_spectra(rydberg_constant, atomic_number, principal_quantum_number_initial,
                       principal_quantum_number_final):
        return (rydberg_constant * atomic_number ** 2) * (
                (1 / principal_quantum_number_final ** 2) - (1 / principal_quantum_number_initial ** 2))

    @staticmethod
    def absorption_spectrum(frequency_emitted, frequency_absorbed, velocity):
        return (frequency_absorbed - frequency_emitted) / frequency_emitted == velocity / speed_of_light

    @staticmethod
    def emission_spectrum(frequency_emitted, frequency_absorbed, velocity):
        return (frequency_absorbed - frequency_emitted) / frequency_absorbed == velocity / speed_of_light

    @staticmethod
    def electron_smath_pin_magnetic_moment(magnetic_moment):
        return magnetic_moment / (Charge.electron * Mass.electron)

    @staticmethod
    def electron_g_factor(g_factor, bohr_magneton):
        return g_factor * bohr_magneton

    @staticmethod
    def atomic_nucleus_g_factor(g_factor, nuclear_magneton):
        return g_factor * nuclear_magneton

    @staticmethod
    def nuclear_decay_constant(half_life):
        return math.log(2) / half_life

    @staticmethod
    def nuclear_decay(initial_quantity, decay_constant, time):
        return initial_quantity * math.exp(-decay_constant * time)
