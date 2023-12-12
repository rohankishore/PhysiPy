from PhysiPy import *


def main():
    print("ran")
    width = 100
    height = 20
    area = equations.area(width, height)
    depth = 40
    answer = depth * area
    print(answer)

if __name__ == '__main__':
    main()
