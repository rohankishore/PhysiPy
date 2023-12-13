from .constants import *


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
    def ideal_gas_law(pressure, v, temperature, gas_constant1):
        return (pressure * v) / (gas_constant1 * temperature)
