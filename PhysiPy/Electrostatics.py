from .constants import Coulombs_constant


class Electrostatics:

    def __init__(self, distance=1, charge=1, charge_1=1, charge_2=1, voltage=1, time=1,
                 current=1, force=1, velocity=1, length=1, charge_enclosed=1,
                 epsilon_0=1, emf=1, electric_field=1, magnetic_field=1,
                 resistance=1, charge_carrier_density=1, thickness=1,
                 cross_sectional_area=1) -> None:
        self.distance = distance
        self.charge = charge
        self.c_1 = charge_1
        self.c_2 = charge_2
        self.voltage = voltage
        self.time = time
        self.current = current
        self.force = force
        self.velocity = velocity
        self.length = length
        self.charge_enclosed = charge_enclosed
        self.epsilon_0 = epsilon_0
        self.emf = emf
        self.electric_field1 = electric_field
        self.magnetic_field1 = magnetic_field
        self.r = resistance
        self.charge_carrier_density = charge_carrier_density
        self.thickness = thickness
        self.cross_sectional_area = cross_sectional_area

    def electric_force(self):
        return Coulombs_constant * self.c_1 * self.c_2 / self.distance ** 2

    def electric_field(self):
        return Coulombs_constant * self.charge / self.distance ** 2

    def electric_potential(self):
        return Coulombs_constant * self.charge / self.distance

    def capacitance(self):
        return self.charge / self.voltage

    def electric_current(self):
        return self.charge / self.time

    def resistance(self):
        return self.voltage / self.current

    def ohms_law(self):
        return self.voltage / self.current

    def coulombs_law(self):
        return self.electric_field1 * self.charge

    def gauss_law(self):
        return self.charge_enclosed / self.epsilon_0

    def faradays_law(self):
        return self.emf * self.time

    def magnetic_field(self):
        return self.force / (self.charge * self.velocity)

    def lorentz_force(self):
        return self.charge * self.velocity * self.magnetic_field1

    def hall_voltage(self):
        return (self.magnetic_field1 * self.current) / \
               (self.charge_carrier_density * self.thickness)

    def drift_velocity(self):
        return self.current / (self.charge * self.charge_carrier_density *
                               self.cross_sectional_area)

    def resistivity(self):
        return self.r * (self.cross_sectional_area / self.length)

    @staticmethod
    def series_resistance(*resistances):
        return sum(resistances)

    @staticmethod
    def parallel_resistance(*resistances):
        return 1 / sum(1 / r for r in resistances)
