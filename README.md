# Thermo Snookered
**By Viraj Patel**
Department of Physics, Imperial College London

## Contents

1. Introduction
2. How to use the main program
3. How to use **testing<span>.py</span>**
4. Important points
5. References

## Introduction


Thermo Snookered is a physics program written in Python 3.8 that uses thermodynamics to simulate gas particles in a circular container. It is restricted to 2 dimensions for simplicity and uses an OOP approach for completeness.

This program contains the following files:
- **collisions_in_circle.py**: Manages collisions
- **data_collection.py**: Collects data in simulation
- **thermo_snookered_v1.py**: All classes for simulation
- **thermo_snookered_frontend.py**: Uses backend to produce graphs
- **testing<span>.py</span>**: Tests conducted during development

The zip file also contains the following image files:
- T-9-1.png: Histogram of distances from centre
- T-9-2.png: Histogram of inter-ball separation
- T-11-1.png: Kinetic energy against time graph
- T-11-2.png: Momentum against time graph
- T-11-3.png: Pressure against temperature graph
- T-12.png: Graph showing how equation of state changes with radius
- T-13.png: Distribution of velocities
- T-14.png: Comparison of equation of state with van der Waals equation


## How to use the main program

1.  First, make sure all of the files stated above are in a single folder on your computer and you have Python 3.8 installed including the following packages:
    - MatPlotLib
    - Numpy
    - Itertools
    - Random

2.  Open all of the Python files except for **testing<span>.py</span>** (this is not required right now). Using **thermo_snookered_frontend.py**, uncomment the code that you would like to run and then execute the file. Please only uncomment one line at a time. The frontend has been carefully written to allow the tasks to be run in separate executions.
3. When running the animation, you will see two sets of axes. The left set of axes has the simulation. On the right is a dynamic graph of dots. The red dot represents the total kinetic energy of the system and should not move. To make sure this was not plotted on the same point in each iteration there are blue dots to represent the individual kinetic energies of the balls that move with each collision. The green dot on the far right represents the total momentum of the system.

**PLEASE NOTE**
All of the procedures in **thermo_snookered_frontend.py** have a *save* parameter, the purpose of this was for me to save the image files. I kept the feature in case saving the plots are required but it is important to note that this will replace the provided image files.

## How to use testing<span>.py</span>

1.  Make sure that the following files are open in your code editor:
    - **thermo_snookered_frontend.py**
    - **thermo_snookered_v1.py**
    - **collisions_in_circle.py**
    - **data_collection.py**

2.  Please follow the instructions in the code. Please uncomment one section at a time to make sure that there are no errors. There are 6 tests that test each class. It is also worth noting that the procedures in the frontend can also be considered as tests for both the functionality of the program and the physics of the system.


## Important points

I used *matplotlib.animation.FuncAnimation* for the animation feature as opposed to the recommended method. This is because it produces a smoother animation and allows me to see if there is any overlapping between collisions or frames.

## References
These are the most important sources I used to help create this program:

[1.](https://scipython.com/blog/two-dimensional-collisions/) christian, *"Two-dimensional Collisions"*, Scipython, 24 June 2019
This source helped me understand and use the *matplotlib.animation.FuncAnimation* function. It also showed me the syntax to use for the linear algebra operations that I needed. 
[2.](https://numpy.org/doc/stable/reference/routines.linalg.html) Numpy, *"Linear algebra (numpy.linalg)"*, 29 June 2020
[3.](https://stackoverflow.com/questions/403421/how-to-sort-a-list-of-objects-based-on-an-attribute-of-the-objects)Triptych, *"How to sort a list of objects based on an attribute of the objects?"*, 31 December 2008
