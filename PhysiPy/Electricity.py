from .constants import Coulombs_constant


class Electricity:
    """
        This class is the Electricity Suit. it contains the functions
        functions for Electricity calculations.
    """
    def __init__(self, voltage=1, current=1, resistance=1,
                 q1=1, q2=1, r=1) -> None:
        self.V = voltage
        self.C = current
        self.R = resistance
        self.k = Coulombs_constant
        self.q1 = q1
        self.q2 = q2
        self.r = r

    def force_electrostatics(self) -> float:
        return (self.k * self.q1 * self.q2) / (self.r ** 2)

    def resistance(self) -> float:
        if self.C != 0 and self.C is not None:
            return self.V / self.C
        else:
            # Handle division by zero
            raise ValueError("Current cannot be zero for resistance calculation.")

    def current(self) -> float:
        if self.R != 0 and self.R is not None:
            return self.V / self.R
        else:
            # Handle division by zero
            raise ValueError("Resistance cannot be zero for current calculation.")

    def voltage(self) -> float:
        return self.C * self.R

    def power(self) -> float:
        return self.V * self.C

    def ohms_law(self):
        if self.R != 0 and self.R is not None:
            return self.V / self.R
        else:
            # Handle division by zero
            raise ValueError("Resistance cannot be zero for ohms law calculation.")

