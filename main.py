from PhysiPy.Electricity import Electricity
# import math


def main():
    print("ran")

    # Setting variables
    voltage = 230
    resistance = 20
    # setting up the class
    e = Electricity(voltage=voltage, resistance=resistance)
    # retrieving and printing the answer
    answer = e.current()

    print(f'{answer:,}')


if __name__ == '__main__':
    main()
