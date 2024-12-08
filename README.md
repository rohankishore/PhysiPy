<h1 align="center">PhysiPy</h1>

<h3 align="center">Physics Equation Solver and Constants for Python</h3>

<div align="center">

  <a href="https://pypi.org/project/PhysiPy-Python/#history">![License](https://img.shields.io/badge/Install-PyPI-blue)</a>
  <a href="https://opensource.org/licenses/MIT">![License](https://img.shields.io/badge/License-MIT-yellow)</a>
  <a href="https://www.fiverr.com/rohancodespy/">![Demo](https://img.shields.io/badge/Fiverr-Hire-green)</a>
</div>

PhysiPy is a powerful and versatile Python library designed to streamline physics calculations and provide easy access to a vast collection of essential physical constants. Whether you are a student, researcher, or an enthusiast seeking to explore the intricacies of the physical world, PhysiPy is an indispensable tool for your scientific endeavors.

With PhysiPy, you can effortlessly perform complex physics computations without the need for extensive manual coding. The library encompasses a wide range of formulas and equations spanning various branches of physics, including mechanics, electromagnetism, thermodynamics, quantum mechanics, and more. From simple kinematic equations to intricate quantum mechanical wavefunctions, PhysiPy has you covered, simplifying the process of implementing these calculations into your code.

One of the key features of PhysiPy is its extensive collection of physical constants. It contains hundreds of well-documented and up-to-date constants that are crucial for numerous calculations and experiments. These constants encompass fundamental values such as the speed of light, Planck's constant, elementary charge, Avogadro's number, and many more. By having this wealth of constants readily available, PhysiPy eliminates the need to search for and manually input these values, ensuring accuracy and efficiency in your calculations.

<br>

<div align="center">
  
PhysiPy API is currently in development by [Vikram Samak](https://github.com/vikramsamak) and can be found [here](https://github.com/vikramsamak/PhsiPy-Api)
</div>

<br>

<h1 align="center">⬇️ Installation</h1>

<div align="center">

```bash
pip install PhysiPy-Python
```

</div>

<h3 align="center">You can also install it via the given .targz file. Here are the steps:</h3>

<br>

- Download the `PhysiPy-1.0.0.tar.gz` file from the `dist` folder in the repository
- Now, copy the path of the downloaded `.tar.gz` file
- Open Terminal and type in the following command:
` pip install <path you've copied>`

<br>

<h1 align="center">⭐ Features in a Glance</h1>

- Over 100 pre-defined Physics Equations (Just substitute the values inside function)
- Over 150+ constants (including Boltzmann Constant, Gravitational Constants, and much more)
- Extremely quick since its uses Numpy

<br>

<h1 align="center">🧑🏻‍💻 Demo Snippets</h1>

<br>

```python
# to calculate resistance

import PhysiPy.Electricity as ec

a = ec.resistance(25, 10)
print(a)

>> 2.5
```

<br>

```python
# to calculate gravitational potential

import PhysiPy.Gravitation as gr

a = gr.gravitational_acceleration(25, 1200)
print(a)

>> 1.1586944444444442e^-15
```

<h1 align="center">❣️ Used by: </h1>

- [Out-of-the-Box-Astronautics-LLC](https://github.com/Out-of-the-Box-Astronautics-LLC/) for their [Lunar Lander project](https://github.com/Out-of-the-Box-Astronautics-LLC/StrongBox)
