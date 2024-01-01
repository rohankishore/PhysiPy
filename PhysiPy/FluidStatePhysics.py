from .constants import pi, speed_of_sound, gas_constant
import math


class FluidStatePhysics:

    def __init__(self, density=1, acceleration_due_to_gravity=1, height=1,
                 surface_tension_coefficient=1, contact_angle=1, radius_of_curvature=1,
                 static_pressure=1, dynamic_pressure=1, gravitational_acceleration=1,
                 flow_rate=1, viscosity=1, length=1, radius=1, surface_tension=1,
                 velocity=1, drag_force=1, particle_radius=1, volume_change=1,
                 initial_volume=1, pressure_change=1, initial_pressure=1,
                 final_pressure=1, final_volume=1, initial_temperature=1,
                 final_temperature=1, initial_moles=1, final_moles=1, pressure=1,
                 v=1, temperature=1):
        self.density = density
        self.acceleration_due_to_gravity = acceleration_due_to_gravity
        self.height = height
        self.surface_tension_coefficient = surface_tension_coefficient
        self.contact_angle = contact_angle
        self.radius_of_curvature = radius_of_curvature
        self.static_pressure = static_pressure
        self.dynamic_pressure = dynamic_pressure
        self.gravitational_acceleration = gravitational_acceleration
        self.flow_rate = flow_rate
        self.viscosity = viscosity
        self.length = length
        self.radius = radius
        self.surface_tension = surface_tension
        self.velocity = velocity
        self.drag_force = drag_force
        self.particle_radius = particle_radius
        self.volume_change = volume_change
        self.initial_volume = initial_volume
        self.pressure_change = pressure_change
        self.initial_pressure = initial_pressure
        self.final_pressure = final_pressure
        self.final_volume = final_volume
        self.initial_temperature = initial_temperature
        self.final_temperature = final_temperature
        self.initial_moles = initial_moles
        self.final_moles = final_moles
        self.pressure = pressure
        self.v = v
        self.temperature = temperature

    def hydrostatic_pressure(self):
        return self.density * self.acceleration_due_to_gravity * self.height

    def surface_tension(self):
        return self.surface_tension_coefficient * math.cos(self.contact_angle) / \
               self.radius_of_curvature

    def capillary_pressure(self):
        return 2 * self.surface_tension / self.radius_of_curvature

    def bernoullis_equation(self):
        return self.static_pressure + self.dynamic_pressure + self.density * \
               self.gravitational_acceleration * self.height

    def poiseuilles_law(self):
        return (pi * self.radius ** 4 * self.flow_rate) / \
               (8 * self.viscosity * self.length)

    def reynolds_number(self):
        return (self.density * self.velocity * self.length) / self.viscosity

    def stokes_law(self):
        return (6 * pi * self.viscosity * self.particle_radius) / self.drag_force

    def mach_number(self):
        return self.velocity / speed_of_sound

    def compressibility_factor(self):
        return self.volume_change / (self.initial_volume * self.pressure_change)

    def boyles_law(self):
        return (self.initial_pressure * self.initial_volume) / \
               (self.final_pressure * self.final_volume)

    def charles_law(self):
        return (self.initial_volume * self.final_temperature) / \
               (self.final_volume * self.initial_temperature)

    def gaylussacs_law(self):
        return (self.initial_pressure * self.final_temperature) / \
               (self.final_pressure * self.initial_temperature)

    def avogadros_law(self):
        return (self.initial_volume * self.final_moles) / \
               (self.final_volume * self.initial_moles)

    def ideal_gas_law(self):
        return (self.pressure * self.v) / (gas_constant * self.temperature)
