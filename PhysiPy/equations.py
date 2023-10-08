import numpy as np
import math


def weight(mass):
    result = mass * g
    return result

def density(mass, volume):
    result = mass/volume
    return result

def volume(L, B, H):
    print("Volume:", L * B * H)

def area(L, B):
    print("Area:", L * B)

class Gravitation:
    def __init__(self):
        super().__init__()

    @staticmethod
    def G(m1, m2, r):
        rsq = r*r
        result = G * ((m1*m2)/rsq)
        return result

    @staticmethod
    def G_Potential(M, r):
        v = (-G * M) / r
        return v

    @staticmethod
    def g_in_depth(depth):
        res = g * (1-(depth/6400))
        return res

    @staticmethod
    def axial_velocity(area_swept, time):
        return area_swept/time

    @staticmethod
    def gravitational_force(mass1, mass2, distance):
        return G * mass1 * mass2 / distance ** 2

    @staticmethod
    def gravitational_potential_energy(mass1, mass2, distance):
        return -(G * mass1 * mass2) / distance

    @staticmethod
    def gravitational_field_strength(mass, distance):
        return G * mass / distance ** 2

    @staticmethod
    def escape_velocity(mass, radius):
        return math.sqrt(2 * G * mass / radius)

    @staticmethod
    def orbital_velocity(mass, radius):
        return math.sqrt(G * mass / radius)

    @staticmethod
    def period_of_orbit(mass, radius):
        return 2 * math.pi * math.sqrt(radius ** 3 / (G * mass))

    @staticmethod
    def gravitational_potential(mass, distance):
        return -(G * mass) / distance

    @staticmethod
    def weight(mass, acceleration_due_to_gravity):
        return mass * acceleration_due_to_gravity

    @staticmethod
    def gravitational_acceleration(mass, distance):
        return G * mass / distance ** 2

    @staticmethod
    def keplers_third_law(period, semi_major_axis):
        return (period ** 2) / (semi_major_axis ** 3)

class Nlm:
    def __init__(self):
        super().__init__()

    @staticmethod
    def force(mass, acceleration):
        return mass * acceleration

    @staticmethod
    def momentum(mass, velocity):
        return mass * velocity

    @staticmethod
    def recoil_velocity(massofbullet, massofgun, initalvelocity):
        return (massofbullet*initalvelocity)/massofgun

class Mechanics:
    def __init__(self):
        super().__init__()

    @staticmethod
    def velocity(initial_velocity, acceleration, time):
        return initial_velocity + acceleration * time

    @staticmethod
    def displacement(initial_velocity, acceleration, time):
        return initial_velocity * time + 0.5 * acceleration * time ** 2

    @staticmethod
    def acceleration(final_velocity, initial_velocity, time):
        return (final_velocity - initial_velocity) / time

    @staticmethod
    def uniform_accelerated_motion(initial_velocity, time, acceleration):
        displacement = initial_velocity * time + 0.5 * acceleration * time ** 2
        final_velocity = initial_velocity + acceleration * time
        return displacement, final_velocity

    @staticmethod
    def force(mass, acceleration):
        return mass * acceleration

    @staticmethod
    def work(force, displacement, angle=0):
        return force * displacement * math.cos(math.radians(angle))

    @staticmethod
    def kinetic_energy(mass, velocity):
        return 0.5 * mass * velocity ** 2

    @staticmethod
    def potential_energy(mass, height, gravitational_field_strength):
        return mass * height * gravitational_field_strength

    @staticmethod
    def power(work, time):
        return work / time

    @staticmethod
    def momentum(mass, velocity):
        return mass * velocity

    @staticmethod
    def impulse(force, time):
        return force * time

    @staticmethod
    def circular_velocity(radius, period):
        return 2 * math.pi * radius / period

    @staticmethod
    def centripetal_acceleration(radius, velocity):
        return velocity ** 2 / radius

    @staticmethod
    def torque(force, lever_arm):
        return force * lever_arm

    @staticmethod
    def angular_velocity(angular_displacement, time):
        return angular_displacement / time

    @staticmethod
    def angular_acceleration(angular_velocity, time):
        return angular_velocity / time

class Errors:
    def __init__(self):
        super().__init__()

    @staticmethod
    def error_muldiv(a, b, c, d):
        x = ((a / c) + (b / d))
        return x

    @staticmethod
    def error_addsub(a, b):
        p = (a + b)
        return p

    @staticmethod
    def percentage_error(M, E):
        print("Percentage Error:", (M / E) * 100, "%")

    @staticmethod
    def absolute_error(a, n):
        if len(a) == n:
            A = np.array(a)
            S = np.sum(A)
            E = S / n
            return E
        else:
            return ValueError

    @staticmethod
    def meanabsolute_error(a, n):
        if len(a) == n:
            b = []
            for i in a:
                b.append(abs(i))
            A = np.array(b)
            S = np.sum(A)
            M = S / n
            return M
        else:
            print("error try again")

class Electricity:
    def __init__(self):
        super().__init__()

    @staticmethod
    def force_electrostatics(q1, q2, r):
        r = r*r
        y = (q1 * q2) / r
        x = permittivity_of_free_space * y
        return x

    @staticmethod
    def resistance(voltage, current):
        p = (voltage/ current)
        return p

    @staticmethod
    def current(voltage, resistance):
        return voltage/resistance

    @staticmethod
    def voltage(current, resistance):
        return current*resistance

    @staticmethod
    def power(voltage, current):
        return voltage * current

class Electrostatics:
    @staticmethod
    def electric_force(charge1, charge2, distance):
        k = 9e9  # Coulomb's constant
        return k * charge1 * charge2 / distance ** 2

    @staticmethod
    def electric_field(charge, distance):
        k = 9e9  # Coulomb's constant
        return k * charge / distance ** 2

    @staticmethod
    def electric_potential(charge, distance):
        k = 9e9  # Coulomb's constant
        return k * charge / distance

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

class Waves:
    @staticmethod
    def wave_velocity(frequency, wavelength):
        return frequency * wavelength

    @staticmethod
    def angular_frequency(frequency):
        return 2 * pi * frequency

    @staticmethod
    def wave_period(frequency):
        return 1 / frequency

    @staticmethod
    def wave_number(wavelength):
        return 2 * pi / wavelength

    @staticmethod
    def wave_speed(wavelength, period):
        return wavelength / period

    @staticmethod
    def longitudinal_wave_speed(frequency, wavelength):
        return frequency * wavelength

    @staticmethod
    def intensity(power, area):
        return power / area

    @staticmethod
    def sound_intensity(sound_power, area):
        return sound_power / area

    @staticmethod
    def sound_level(intensity):
        return 10 * math.log10(intensity / (10 ** -12))

    @staticmethod
    def beats_frequency(frequency1, frequency2):
        return abs(frequency1 - frequency2)

    @staticmethod
    def beats_period(frequency1, frequency2):
        return 1 / abs(frequency1 - frequency2)

    @staticmethod
    def doppler_effect(observed_frequency, source_frequency, wave_velocity, observer_velocity, source_velocity):
        return observed_frequency / source_frequency == (
                wave_velocity + observer_velocity) / (wave_velocity + source_velocity)

    @staticmethod
    def refractive_index(real_depth, apparent_depth):
        return real_depth / apparent_depth

class Electromagnetism:
    @staticmethod
    def voltage(current, resistance):
        return current * resistance

    @staticmethod
    def power(current, voltage):
        return current * voltage

    @staticmethod
    def resistance(voltage, current):
        return voltage / current

    @staticmethod
    def resistivity(resistance, area, length):
        return (resistance * area) / length

    @staticmethod
    def electric_field(voltage, distance):
        return voltage / distance

    @staticmethod
    def electric_potential_energy(charge, voltage):
        return charge * voltage

    @staticmethod
    def electric_power(current, voltage):
        return current * voltage

    @staticmethod
    def magnetic_field_strength(force, length, current):
        return force / (length * current)

    @staticmethod
    def magnetic_flux_density(magnetic_flux, area):
        return magnetic_flux / area

    @staticmethod
    def magnetic_flux(magnetic_flux_density, area):
        return magnetic_flux_density * area

    @staticmethod
    def magnetic_force(magnetic_field_strength, length, current):
        return magnetic_field_strength * length * current

    @staticmethod
    def lorentz_force(charge, magnetic_field_strength, velocity):
        return charge * magnetic_field_strength * velocity

    @staticmethod
    def hall_voltage(magnetic_field_strength, current, charge_density, electron_density):
        return (magnetic_field_strength * current) / (charge_density * electron_density)

    @staticmethod
    def faradays_law(induced_emf, time, number_of_turns):
        return induced_emf * time / number_of_turns

    @staticmethod
    def self_inductance(inductance, current_rate_of_change):
        return inductance * current_rate_of_change

    @staticmethod
    def mutual_inductance(induced_emf, current_rate_of_change):
        return induced_emf / current_rate_of_change

    @staticmethod
    def coulombs_law(charge1, charge2, distance):
        return (1 / (4 * pi * epsilon0)) * ((charge1 * charge2) / distance ** 2)

    @staticmethod
    def capacitance(charge, voltage):
        return charge / voltage

    @staticmethod
    def electric_field_strength(charge, distance):
        return (1 / (4 * pi * epsilon0)) * (charge / distance ** 2)

    @staticmethod
    def electric_flux_density(charge, area):
        return charge / area

    @staticmethod
    def ohms_law(voltage, resistance):
        return voltage / resistance

    @staticmethod
    def magnetic_flux_quantum(magnetic_flux):
        return magnetic_flux / (vonklitzing_constant / (2 * pi))

    @staticmethod
    def de_broglie_wavelength(momentum):
        return plancks_constant / momentum

    @staticmethod
    def compton_wavelength_change(initial_wavelength, scattering_angle):
        return initial_wavelength - (compton_wavelength * (1 - math.cos(scattering_angle)))

    @staticmethod
    def bohr_orbit_radius(principal_quantum_number):
        return (finestructure_constant ** 2) * (permittivity_of_free_space * (plancks_constant ** 2)) / (
                (math.pi * elementary_charge ** 2) * (principal_quantum_number ** 2) * (speed_of_light ** 2))

    @staticmethod
    def fermi_energy(fermi_level):
        return fermi_level * elementary_charge

    @staticmethod
    def de_broglie_wavelength_matter_wave(mass, velocity):
        return plancks_constant / (mass * velocity)

    @staticmethod
    def number_of_ions(charge, elementary_charge):
        return charge / elementary_charge

    @staticmethod
    def acceleration_due_to_gravity(gravitational_field_strength, mass):
        return gravitational_field_strength * mass

    @staticmethod
    def magnetic_field_inside_a_solenoid(mu0, number_of_turns, current, length):
        return (mu0 * number_of_turns * current) / length

    @staticmethod
    def magnetic_field_around_a_wire(mu0, current, radius):
        return (mu0 * current) / (2 * pi * radius)

    @staticmethod
    def induced_emf(change_in_magnetic_flux, time):
        return -change_in_magnetic_flux / time

    @staticmethod
    def maxwell_equation(magnetic_field, magnetic_flux_density):
        return np.cross(magnetic_field, magnetic_flux_density) - permittivity_of_free_space * magnetic_field

    @staticmethod
    def snells_law(n1, n2, angle_of_incidence):
        return (n1 * math.sin(angle_of_incidence)) / n2

    @staticmethod
    def critical_angle(n1, n2):
        return math.asin(n2 / n1)

    @staticmethod
    def photoelectric_effect_work_function(hv, threshold_frequency):
        return hv - threshold_frequency

    @staticmethod
    def photoelectric_effect_max_velocity(work_function, stopping_potential):
        return math.sqrt((2 * elementary_charge * (stopping_potential - work_function)) / Subatomic.Mass.electron)

    @staticmethod
    def decay_law(initial_quantity, decay_constant, time):
        return initial_quantity * math.exp(-decay_constant * time)

    @staticmethod
    def half_life(decay_constant):
        return math.log(2) / decay_constant

    @staticmethod
    def nuclear_binding_energy(mass_defect):
        return mass_defect * speed_of_light ** 2


    @staticmethod
    def radioactive_decay(initial_quantity, final_quantity, time):
        return -initial_quantity * math.log(final_quantity / initial_quantity) / time

    @staticmethod
    def einstein_mass_energy_equivalence(mass):
        return mass * speed_of_light ** 2

    @staticmethod
    def coulomb_law(charge1, charge2, distance):
        return (1 / (4 * pi * epsilon0)) * ((charge1 * charge2) / distance ** 2)

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
    def first_law_thermodynamics(heat_added, work_done, change_in_internal_energy):
        return heat_added + work_done == change_in_internal_energy

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
    def average_kinetic_energy(molar_mass, temperature):
        return (3 / 2) * boltzmann_constant * temperature

    @staticmethod
    def speed_of_sound(molar_mass, temperature, gamma):
        return math.sqrt(gamma * boltzmann_constant * temperature / molar_mass)

    @staticmethod
    def specific_heat_capacity(mass, temperature_change, heat_transfer):
        return heat_transfer / (mass * temperature_change)

    @staticmethod
    def latent_heat(mass, heat_transfer):
        return heat_transfer / mass

class QuantumMechanics:
    @staticmethod
    def uncertainty_principle(position_uncertainty, momentum_uncertainty):
        return position_uncertainty * momentum_uncertainty >= plancks_constant / (4 * pi)

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
                speed_of_light ** 2 * (math.exp((plancks_constant * frequency) / (boltzmann_constant * temperature)) - 1))

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
        return energy * lifetime >= plancks_constant / (4 * pi)

    @staticmethod
    def heisenberg_uncertainty_position_momentum(position_uncertainty, momentum_uncertainty):
        return position_uncertainty * momentum_uncertainty >= plancks_constant / 2

    @staticmethod
    def atomic_orbital_radius(principal_quantum_number):
        return (4 * pi * permittivity_of_free_space * (plancks_constant ** 2) * (
                principal_quantum_number ** 2)) / (Subatomic.Mass.electron * (elementary_charge ** 2))


    @staticmethod
    def fine_structure_constant():
        return (elementary_charge ** 2) / (4 * pi * permittivity_of_free_space * plancks_constant * speed_of_light)

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
                (plancks_constant ** 3) * (speed_of_light ** 2) * (math.exp(energy / (boltzmann_constant * temperature)) - 1))

    @staticmethod
    def plank_law_intensity(frequency, temperature):
        return (2 * plancks_constant * frequency ** 3) / (
                speed_of_light ** 2 * (math.exp((plancks_constant * frequency) / (boltzmann_constant * temperature)) - 1))

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
    def electron_spin_magnetic_moment(magnetic_moment):
        return magnetic_moment / (Subatomic.Charge.electron * Subatomic.Mass.electron)

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
        return (electron_density * (plancks_constant ** 2)) / (2 * Subatomic.Mass.electron) + work_function

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
    def intrinsic_carrier_concentration(donor_density, acceptor_density, fermi_energy, intrinsic_fermi_level, temperature):
        return ((donor_density * acceptor_density) ** 0.5) * math.exp(
            -(fermi_energy - intrinsic_fermi_level) / (2 * boltzmann_constant * temperature))

    @staticmethod
    def depletion_layer_width(fermi_energy, potential_difference, electron_density):
        return math.sqrt((2 * permittivity_of_free_space * potential_difference) / (
                Subatomic.Charge.electron * electron_density)) / (2 * fermi_energy)

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
    def thermal_conductivity(heat_flux, temperature_gradient, thickness):
        return heat_flux * thickness / (temperature_gradient * area)

    @staticmethod
    def superconductivity_critical_temperature(critical_field):
        return (critical_field ** 2) * (permittivity_of_free_space / (2 * Subatomic.Mass.electron * Subatomic.Charge.electron))

    @staticmethod
    def superconductivity_coherence_length(hamiltonian_operator, fermi_energy):
        return (plancks_constant / (2 * math.pi)) * (
                (3 * permittivity_of_free_space * (hamiltonian_operator ** 2)) / (2 * Subatomic.Mass.electron * fermi_energy)) ** 0.5

class FluidStatePhysics:
    @staticmethod
    def hydrostatic_pressure(density, acceleration_due_to_gravity, height):
        return density * acceleration_due_to_gravity * height

    @staticmethod
    def surface_tension(surface_tension_coefficient, contact_angle, radius_of_curvature):
        return surface_tension_coefficient * math.cos(contact_angle) / radius_of_curvature

    @staticmethod
    def capillary_pressure(surface_tension, radius_of_curvature):
        return 2 * surface_tension / radius_of_curvature

    @staticmethod
    def bernoullis_equation(static_pressure, dynamic_pressure, height, density, gravitational_acceleration):
        return static_pressure + dynamic_pressure + density * gravitational_acceleration * height

    @staticmethod
    def poiseuilles_law(flow_rate, viscosity, length, radius):
        return (math.pi * radius ** 4 * flow_rate) / (8 * viscosity * length)

    @staticmethod
    def reynolds_number(density, velocity, length, viscosity):
        return (density * velocity * length) / viscosity

    @staticmethod
    def stokes_law(drag_force, viscosity, particle_radius):
        return (6 * math.pi * viscosity * particle_radius) / drag_force

    @staticmethod
    def mach_number(speed_of_sound, velocity):
        return velocity / speed_of_sound

    @staticmethod
    def compressibility_factor(volume_change, initial_volume, pressure_change):
        return volume_change / (initial_volume * pressure_change)

    @staticmethod
    def boyles_law(initial_pressure, initial_volume, final_pressure, final_volume):
        return (initial_pressure * initial_volume) / (final_pressure * final_volume)

    @staticmethod
    def charles_law(initial_volume, final_volume, initial_temperature, final_temperature):
        return (initial_volume * final_temperature) / (final_volume * initial_temperature)

    @staticmethod
    def gaylussacs_law(initial_pressure, final_pressure, initial_temperature, final_temperature):
        return (initial_pressure * final_temperature) / (final_pressure * initial_temperature)

    @staticmethod
    def avogadros_law(initial_volume, final_volume, initial_moles, final_moles):
        return (initial_volume * final_moles) / (final_volume * initial_moles)

    @staticmethod
    def ideal_gas_law(pressure, volume, temperature, gas_constant):
        return (pressure * volume) / (gas_constant * temperature)