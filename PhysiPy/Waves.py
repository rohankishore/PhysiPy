from .constants import pi
import math


class Waves:

    def __init__(self, frequency=1, wavelength=1, period=1, power=1, area=1,
                 sound_power=1, intensity=1, frequency1=1, frequency2=1,
                 velocity_wave=1, velocity_observer=1, velocity_source=1,
                 real_depth=1, apparent_depth=1):
        self.frequency = frequency
        self.wavelength = wavelength
        self.period = period
        self.power = power
        self.area = area
        self.sound_power = sound_power
        self.intensity = intensity
        self.frequency1 = frequency1
        self.frequency2 = frequency2
        self.velocity_wave = velocity_wave
        self.velocity_observer = velocity_observer
        self.velocity_source = velocity_source
        self.real_depth = real_depth
        self.apparent_depth = apparent_depth

    def wave_velocity(self):
        return self.frequency * self.wavelength

    def angular_frequency(self):
        return 2 * pi * self.frequency

    def wave_period(self):
        return 1 / self.frequency

    def wave_number(self):
        return 2 * pi / self.wavelength

    def wave_speed(self):
        return self.wavelength / self.period

    def longitudinal_wave_speed(self):
        return self.frequency * self.wavelength

    def intensity(self):
        return self.power / self.area

    def sound_intensity(self):
        return self.sound_power / self.area

    def sound_level(self):
        return 10 * math.log10(self.intensity / (10 ** -12))
    
    def beats_frequency(self):
        return abs(self.frequency1 - self.frequency2)

    def beats_period(self):
        return 1 / abs(self.frequency1 - self.frequency2)

    def doppler_frequency(self):
        return self.frequency * (self.velocity_wave + self.velocity_observer) / \
               (self.velocity_wave - self.velocity_source)

    def doppler_wavelength(self):
        return self.wavelength * (self.velocity_wave - self.velocity_source) / \
               (self.velocity_wave + self.velocity_observer)

    def refractive_index(self):
        return self.real_depth / self.apparent_depth


"""
You have both an instance variable named intensity and a method with the same name (intensity). 
This can lead to confusion and potential bugs. Consider renaming either the instance variable or 
the method to avoid conflicts.

python

def wave_intensity(self):
    return self.power / self.area

def sound_intensity(self):
    return self.sound_power / self.area

Redundant Intensity Method:
You have two methods calculating intensity: wave_intensity and sound_intensity. Since these seem to serve 
similar purposes, you might want to consider consolidating them into a single method or providing a clearer 
distinction between the two.

def intensity(self, sound=False):
    if sound:
        return self.sound_power / self.area
    else:
        return self.power / self.area
"""
