from .constants import pi, epsilon0, vonklitzing_constant, plancks_constant, \
    compton_wavelength
from .constants import finestructure_constant, Coulombs_constant, elementary_charge, \
    speed_of_light
import math
import numpy as np
from .Mass import Mass


class Electromagnetism:

    def __init__(self, current=1, voltage=1, resistance=1, area=1, length=1, distance=1,
                 charge=1, force=1, magnetic_flux=1, initial_wavelength=1,
                 scattering_angle=1, magnetic_flux_density=1, magnetic_field_strength=1,
                 velocity=1, charge_density=1, electron_density=1, induced_emf=1,
                 time=1, number_of_turns=1, inductance=1, current_rate_of_change=1,
                 charge1=1, charge2=1, principal_quantum_number=1, momentum=1,
                 fermi_level=1, mass=1, elementary_charge1=1, change_in_magnetic_flux=1,
                 magnetic_field=1, gravitational_field_strength=1, mu0=1, n1=1, n2=1,
                 angle_of_incidence=1, hv=1, threshold_frequency=1, work_function=1,
                 stopping_potential=1, initial_quantity=1, decay_constant=1,
                 mass_defect=1, final_quantity=1) -> None:
        self.current = current
        self.voltage = voltage
        self.resistance = resistance
        self.area = area
        self.length = length
        self.distance = distance
        self.charge = charge
        self.force = force
        self.magnetic_flux = magnetic_flux
        self.initial_wavelength = initial_wavelength
        self.scattering_angle = scattering_angle
        self.magnetic_flux_density = magnetic_flux_density
        self.magnetic_field_strength = magnetic_field_strength
        self.velocity = velocity
        self.charge_density = charge_density
        self.electron_density = electron_density
        self.induced_emf = induced_emf
        self.time = time
        self.number_of_turns = number_of_turns
        self.inductance = inductance
        self.current_rate_of_change = current_rate_of_change
        self.charge1 = charge1
        self.charge2 = charge2
        self.principal_quantum_number = principal_quantum_number
        self.momentum = momentum
        self.fermi_level = fermi_level
        self.mass = mass
        self.elementary_charge1 = elementary_charge1
        self.change_in_magnetic_flux = change_in_magnetic_flux
        self.magnetic_field = magnetic_field
        self.gravitational_field_strength = gravitational_field_strength
        self.mu0 = mu0
        self.n1 = n1
        self.n2 = n2
        self.angle_of_incidence = angle_of_incidence
        self.hv = hv
        self.threshold_frequency = threshold_frequency
        self.work_function = work_function
        self.stopping_potential = stopping_potential
        self.initial_quantity = initial_quantity
        self.decay_constant = decay_constant
        self.mass_defect = mass_defect
        self.final_quantity = final_quantity

    def resistivity(self):
        return (self.resistance * self.area) / self.length

    def electric_field(self):
        return self.voltage / self.distance

    def electric_potential_energy(self):
        return self.charge * self.voltage

    def magnetic_field_strength(self):
        return self.force / (self.length * self.current)

    def magnetic_flux_density(self):
        return self.magnetic_flux / self.area

    def magnetic_flux(self):
        return self.magnetic_flux_density * self.area

    def magnetic_force(self):
        return self.magnetic_field_strength * self.length * self.current

    def lorentz_force(self):
        return self.charge * self.magnetic_field_strength * self.velocity

    def hall_voltage(self):
        return (self.magnetic_field_strength * self.current) / (self.charge_density *
                                                                self.electron_density)

    def faradays_law(self):
        return self.induced_emf * self.time / self.number_of_turns

    def self_inductance(self):
        return self.inductance * self.current_rate_of_change

    def mutual_inductance(self):
        return self.induced_emf / self.current_rate_of_change

    def coulombs_law(self):
        return (1 / (4 * pi * epsilon0)) * ((self.charge1 * self.charge2) /
                                            self.distance ** 2)

    def capacitance(self):
        return self.charge / self.voltage

    def electric_field_strength(self):
        return (1 / (4 * pi * epsilon0)) * (self.charge / self.distance ** 2)

    def electric_flux_density(self):
        return self.charge / self.area

    def magnetic_flux_quantum(self):
        return self.magnetic_flux / (vonklitzing_constant / (2 * pi))

    def de_broglie_wavelength(self):
        return plancks_constant / self.momentum

    def compton_wavelength_change(self):
        return self.initial_wavelength - \
               (compton_wavelength * (1 - math.cos(self.scattering_angle)))

    def bohr_orbit_radius(self):
        return (finestructure_constant ** 2) * (Coulombs_constant *
                (plancks_constant ** 2)) / ((pi * elementary_charge ** 2) *
                (self.principal_quantum_number ** 2) * (speed_of_light ** 2))

    def fermi_energy(self):
        return self.fermi_level * elementary_charge

    def de_broglie_wavelength_matter_wave(self):
        return plancks_constant / (self.mass * self.velocity )

    def number_of_ions(self):
        return self.charge / self.elementary_charge1

    def acceleration_due_to_gravity(self):
        return self.gravitational_field_strength * self.mass

    def magnetic_field_inside_a_solenoid(self):
        return (self.mu0 * self.number_of_turns * self.current) / self.length

    def magnetic_field_around_a_wire(self, radius):
        return (self.mu0 * self.current) / (2 * pi * radius)

    def induced_emf(self):
        return -self.change_in_magnetic_flux / self.time

    def maxwell_equation(self):
        return np.cross(self.magnetic_field, self.magnetic_flux_density) - \
               Coulombs_constant * self.magnetic_field

    def snells_law(self):
        return (self.n1 * math.sin(self.angle_of_incidence)) / self.n2

    def critical_angle(self):
        return math.asin(self.n2 / self.n1)

    def photoelectric_effect_work_function(self):
        return self.hv - self.threshold_frequency

    def photoelectric_effect_max_velocity(self):
        return math.sqrt((2 * elementary_charge *
                        (self.stopping_potential - self.work_function)) / Mass.electron)

    def decay_law(self):
        return self.initial_quantity * math.exp(-self.decay_constant * self.time)

    def half_life(self):
        return math.log(2) / self.decay_constant

    def nuclear_binding_energy(self):
        return self.mass_defect * speed_of_light ** 2

    def radioactive_decay(self):
        return -self.initial_quantity * math.log(self.final_quantity /
                                                 self.initial_quantity) / self.time

    def einstein_mass_energy_equivalence(self):
        return self.mass * speed_of_light ** 2

    def coulomb_law(self):
        return (1 / (4 * pi * epsilon0)) * ((self.charge1 * self.charge2) /
                                            self.distance ** 2)
