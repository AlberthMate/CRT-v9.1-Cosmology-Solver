# CRT-v9.1-Cosmology-Solver
I'm sharing my latest proposal on Dark Energy based on LQG. It's called CRT v9.1 and has just been published on Zenodo. We've achieved results consistent with DESI DR2 by eliminating previous numerical strains. If you're interested in open-source cosmology, check it out here: https://doi.org/10.5281/zenodo.19798497

# CRT-v9.1: Toroidal Recoherence Cycle Model

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19798497.svg)](https://doi.org/10.5281/zenodo.19798497)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

The **Toroidal Recoherence Cycle (CRT) v9.1** is a phenomenological framework for Dark Energy inspired by **Loop Quantum Gravity (LQG)** and Lindblad dynamics. This model explores the hypothesis that the accelerated expansion of the universe is driven by a quantum-to-classical transition (recoherence) of the cosmic volume entropy.

### Key Features of v9.1:
* **Dimensional Reformulation:** Resolves the arithmetic overflow issues ($L/\ell_p^2 \sim 10^{96}$) found in previous versions by implementing a dimensionless system.
* **Saturation Mechanism:** Introduces a mobility-dependent "freeze-in" mechanism that stabilizes the equation of state $w$ toward $-1$.
* **Observational Consistency:** Provides quantitative results for $w_0$ and $\rho_\Lambda$ compatible with **DESI DR2 (2024/2026)** and **Planck 2018** data within $1.6\sigma$.

---

## Getting Started

This repository contains the core RK45 numerical integrator used to solve the coupled cosmological equations of the CRT model and generate the evolution of the equation of state $w(z)$.

### Prerequisites

To run the simulations, you will need **Python 3.8+** and the following scientific libraries:

* **NumPy:** For high-performance numerical array operations.
* **SciPy:** For the adaptive Runge-Kutta (RK45) integration methods.
* **Matplotlib:** For generating the cosmological evolution plots ($w(z)$ and $\Lambda(z)$).

### Installation

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/your-username/CRT-v9.1.git](https://github.com/your-username/CRT-v9.1.git)
    cd CRT-v9.1
    ```

2.  **Install dependencies:**
    Using `pip`:
    ```bash
    pip install numpy scipy matplotlib
    ```

---

