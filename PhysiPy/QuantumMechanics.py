from .constants import Coulombs_constant, plancks_constant, speed_of_light
from .constants import elementary_charge, pi, boltzmann_constant
import math
from .Mass import Mass
from .Charge import Charge


class QuantumMechanics:

    def __init__(self, frequency=1, temperature=1, wave_function=1,
                 position_uncertainty=1, momentum_uncertainty=1,
                 hamiltonian_operator=1, energy=1, energy_barrier=1, particle_mass=1,
                 particle_energy=1, fermi_level=1, chemical_potential=1,
                 wavelength=1, momentum=1, lifetime=1, principal_quantum_number=1,
                 number_of_quanta=1, time=1, rydberg_constant=1, atomic_number=1,
                 principal_quantum_number_initial=1, principal_quantum_number_final=1,
                 frequency_emitted=1, frequency_absorbed=1, velocity=1,
                 magnetic_moment=1, g_factor=1, bohr_magneton=1, nuclear_magneton=1,
                 half_life=1, initial_quantity=1, decay_constant=1):
        self.frequency = frequency
        self.temperature = temperature
        self.wave_function = wave_function
        self.position_uncertainty = position_uncertainty
        self.momentum_uncertainty = momentum_uncertainty
        self.hamiltonian_operator = hamiltonian_operator
        self.energy = energy
        self.energy_barrier = energy_barrier
        self.particle_mass = particle_mass
        self.particle_energy = particle_energy
        self.fermi_level = fermi_level
        self.chemical_potential = chemical_potential
        self.wavelength = wavelength
        self.momentum = momentum
        self.lifetime = lifetime
        self.principal_quantum_number = principal_quantum_number
        self.number_of_quanta = number_of_quanta
        self.time = time
        self.rydberg_constant= rydberg_constant
        self.atomic_number = atomic_number
        self.principal_quantum_number_initial = principal_quantum_number_initial
        self.principal_quantum_number_final = principal_quantum_number_final
        self.frequency_emitted = frequency_emitted
        self.frequency_absorbed = frequency_absorbed
        self.velocity = velocity
        self.magnetic_moment = magnetic_moment
        self.g_factor = g_factor
        self.bohr_magneton = bohr_magneton
        self.nuclear_magneton = nuclear_magneton
        self.half_life = half_life
        self.initial_quantity = initial_quantity
        self.decay_constant = decay_constant

    def uncertainty_principle(self):
        return self.position_uncertainty * self.momentum_uncertainty >= \
               plancks_constant / (4 * pi)

    def satisfy_schrodingers_equation(self) -> bool:
        return self.hamiltonian_operator * self.wave_function == self.energy * \
               self.wave_function

    def probability_density(self):
        return (self.wave_function.conjugate() * self.wave_function).real

    def tunneling_probability(self):
        return math.exp((-2 * math.sqrt(2 * self.particle_mass * self.energy_barrier) /
                         plancks_constant) * math.sqrt(self.particle_energy -
                                                       self.energy_barrier))

    def black_body_radiation_intensity(self):
        return (2 * plancks_constant * self.frequency ** 3) / (
                speed_of_light ** 2 * (math.exp((plancks_constant * self.frequency) /
                                        (boltzmann_constant * self.temperature)) - 1))

    def black_body_radiation_power(self):
        return (2 * pi ** 5 * boltzmann_constant ** 4 * self.temperature ** 4) / \
               (15 * plancks_constant ** 3 * speed_of_light ** 2)

    def fermi_dirac_distribution(self):
        return 1 / (math.exp((self.energy - self.fermi_level) / (boltzmann_constant *
                                                                 self.temperature)) + 1)

    def bose_einstein_distribution(self):
        return 1 / (math.exp((self.energy - self.chemical_potential) /
                             (boltzmann_constant * self.temperature)) - 1)

    def exhibits_wave_particle_duality(self) -> bool:
        return self.wavelength * self.momentum == plancks_constant

    def heisenberg_uncertainty_energy_lifetime(self):
        return self.energy * self.lifetime >= plancks_constant / (4 * pi)

    def satisfies_heisenberg_uncertainty_position_momentum(self) -> bool:
        return self.position_uncertainty * self.momentum_uncertainty >= \
               plancks_constant / 2

    def atomic_orbital_radius(self):
        return (4 * pi * Coulombs_constant * (plancks_constant ** 2) * (
                self.principal_quantum_number ** 2)) / (Mass.electron *
                                                        (elementary_charge ** 2))

    def fine_structure_constant(self):
        return (elementary_charge ** 2) / (4 * pi * Coulombs_constant *
                                           plancks_constant * speed_of_light)

    def de_broglie_wavelength_photon(self):
        return speed_of_light / self.frequency

    def de_broglie_wavelength_particle(self):
        return plancks_constant / self.momentum

    def plank_distribution_radiation(self):
        return (2 * plancks_constant * (self.energy ** 3)) / (speed_of_light ** 2 *
                (math.exp((plancks_constant * self.energy) / (boltzmann_constant *
                                                              self.temperature)) - 1))

    def plank_distribution(self):
        return (2 * self.energy ** 2) / (
                (plancks_constant ** 3) * (speed_of_light ** 2) *
                (math.exp(self.energy / (boltzmann_constant * self.temperature)) - 1))

    def plank_law_intensity(self):
        return (2 * plancks_constant * self.frequency ** 3) / (
                speed_of_light ** 2 * (math.exp((plancks_constant * self.frequency) /
                                        (boltzmann_constant * self.temperature)) - 1))

    def plank_law_power(self):
        return (2 * pi * (boltzmann_constant ** 4) * self.temperature ** 4) / \
               (15 * (plancks_constant ** 3) * (speed_of_light ** 2))

    def follows_einstein_light_quanta(self) -> bool:
        return self.energy == plancks_constant * self.frequency

    def einstein_light_intensity(self):
        return self.number_of_quanta / self.time

    def atomic_spectra(self):
        return (self.rydberg_constant * self.atomic_number ** 2) * (
                (1 / self.principal_quantum_number_final ** 2) -
                (1 / self.principal_quantum_number_initial ** 2))

    def has_absorption_spectrum(self) -> bool:
        return (self.frequency_absorbed - self.frequency_emitted) / \
               self.frequency_emitted == self.velocity / speed_of_light

    def has_emission_spectrum(self) -> bool:
        return (self.frequency_absorbed - self.frequency_emitted) / \
               self.frequency_absorbed == self.velocity / speed_of_light

    def electron_smath_pin_magnetic_moment(self):
        return self.magnetic_moment / (Charge.electron * Mass.electron)

    def electron_g_factor(self):
        return self.g_factor * self.bohr_magneton

    def atomic_nucleus_g_factor(self):
        return self.g_factor * self.nuclear_magneton

    def nuclear_decay_constant(self):
        return math.log(2) / self.half_life

    def nuclear_decay(self):
        return self.initial_quantity * math.exp(-self.decay_constant * self.time)
