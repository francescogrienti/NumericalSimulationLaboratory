# Numerical Simulation Laboratory (NSL)

This repository contains the projects and exercises developed for the **Numerical Simulation Laboratory** course at the University of Milan (Università degli Studi di Milano). The course focuses on the application of computational methods and simulation techniques to physical systems.

## Overview

The project is organized into 12 lectures, each covering fundamental topics in numerical physics, ranging from basic statistical methods to advanced machine learning applications. Most simulations are implemented in **C++**, while data analysis, visualization, and deep learning tasks are performed using **Python** and **Jupyter Notebooks**.

## Project Structure

The repository is structured by lecture, with each directory containing the source code, data analysis notebooks, and results.

- **[LSN_lecture_01](./LSN_lecture_01/)**: Random Number Generation & Statistics
    - Testing Pseudo-Random Number Generators (PRNG).
    - Central Limit Theorem validation.
    - Buffon’s Needle experiment for estimating $\pi$.
- **[LSN_lecture_02](./LSN_lecture_02/)**: Monte Carlo Integration & Random Walks
    - Monte Carlo integration using Importance Sampling.
    - Simulation of 3D Random Walks in discrete and continuous space.
- **[LSN_lecture_03](./LSN_lecture_03/)**: Computational Finance
    - European Option Pricing (Call/Put) using Black-Scholes theory.
    - Comparison between direct sampling and discrete path sampling.
- **[LSN_lecture_04](./LSN_lecture_04/)**: Molecular Dynamics (MD)
    - Simulation of Lennard-Jones systems (Argon) in Solid, Liquid, and Gas phases.
    - Implementation of the Verlet algorithm and system equilibration.
- **[LSN_lecture_05](./LSN_lecture_05/)**: Variational Monte Carlo
    - Metropolis algorithm for sampling the Hydrogen atom wavefunctions ($1s$ and $2p$).
- **[LSN_lecture_06](./LSN_lecture_06/)**: Ising Model
    - 1D Ising Model simulation using Metropolis and Gibbs sampling.
    - Phase transitions and thermodynamic properties (Internal Energy, Heat Capacity, Susceptibility).
- **[LSN_lecture_07](./LSN_lecture_07/)**: Advanced MD/MC & Radial Distribution Function
    - NVT Ensemble simulations.
    - Calculation of the Radial Distribution Function $g(r)$ and its comparison across phases.
- **[LSN_lecture_08](./LSN_lecture_08/)**: Quantum Monte Carlo
    - Variational Monte Carlo for the ground state of a 1D quantum system.
    - Introduction to Path Integral Ground State (PIGS) and Path Integral Monte Carlo (PIMC).
- **[LSN_lecture_09](./LSN_lecture_09/)**: Genetic Algorithms (GA)
    - Solving the Traveling Salesman Problem (TSP) using Genetic Algorithms.
- **[LSN_lecture_10](./LSN_lecture_10/)**: Parallel Computing
    - Scaling the TSP solver using MPI (Message Passing Interface) with a "Multiple Islands" approach.
- **[LSN_lecture_11](./LSN_lecture_11/)**: Machine Learning
    - Supervised learning applications: Linear and Polynomial regression using Neural Networks.
- **[LSN_lecture_12](./LSN_lecture_12/)**: Deep Learning
    - Handwritten digit recognition (MNIST dataset) using Convolutional Neural Networks (CNNs) in Keras/TensorFlow.

## Technologies Used

- **Programming Languages**: C++, Python
- **Libraries & Frameworks**:
    - **C++**: MPI (for parallelization), standard mathematical libraries.
    - **Python**: NumPy, SciPy, Matplotlib, Pandas, Keras, TensorFlow.
- **Tools**: Makefile, Jupyter Notebook.

## Getting Started

### Prerequisites

To compile the C++ code, you will need a modern C++ compiler (e.g., `g++`) and optionally `openmpi` for the parallel exercises. For the analysis, a Python environment with the following packages is recommended:
```bash
pip install numpy matplotlib tensorflow jupyter
```

### Running Simulations

Each lecture directory typically contains a `Makefile`. You can compile the simulation using:
```bash
make
./main.exe
```
After running the simulation, the results (usually `.dat` files) can be analyzed using the provided Jupyter Notebooks (`.ipynb`).

## Author
**Francesco Grienti** - *University of Milan*

## License
This project was developed for educational purposes as part of the Numerical Simulation Laboratory course.
