from .constants import permittivity_of_free_space


class Electricity:

    @staticmethod
    def force_electrostatics(q1, q2, r):
        r = r * r
        y = (q1 * q2) / r
        x = permittivity_of_free_space * y
        return x

    @staticmethod
    def resistance(voltage, current):
        p = (voltage / current)
        return p

    @staticmethod
    def current(voltage, resistance):
        return voltage / resistance

    @staticmethod
    def voltage(current, resistance):
        return current * resistance

    @staticmethod
    def power(voltage, current):
        return voltage * current
