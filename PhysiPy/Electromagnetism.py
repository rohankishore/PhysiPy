from .constants import *
import math
import numpy as np
from .Subatomic import *


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
        return (1 / (4 * math.pi * epsilon0)) * ((charge1 * charge2) / distance ** 2)

    @staticmethod
    def capacitance(charge, voltage):
        return charge / voltage

    @staticmethod
    def electric_field_strength(charge, distance):
        return (1 / (4 * math.pi * epsilon0)) * (charge / distance ** 2)

    @staticmethod
    def electric_flux_density(charge, area):
        return charge / area

    @staticmethod
    def ohms_law(voltage, resistance):
        return voltage / resistance

    @staticmethod
    def magnetic_flux_quantum(magnetic_flux):
        return magnetic_flux / (vonklitzing_constant / (2 * math.pi))

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
        return (mu0 * current) / (2 * math.pi * radius)

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
        return math.sqrt((2 * elementary_charge * (stopping_potential - work_function)) / SubatomiMass.electron)

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
        return (1 / (4 * math.pi * epsilon0)) * ((charge1 * charge2) / distance ** 2)
