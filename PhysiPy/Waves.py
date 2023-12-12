from .constants import *
import math


class Waves:
    @staticmethod
    def wave_velocity(frequency, wavelength):
        return frequency * wavelength

    @staticmethod
    def angular_frequency(frequency):
        return 2 * math.pi * frequency

    @staticmethod
    def wave_period(frequency):
        return 1 / frequency

    @staticmethod
    def wave_number(wavelength):
        return 2 * math.pi / wavelength

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
