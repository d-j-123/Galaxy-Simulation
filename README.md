# Galaxy Simulation

A 2D gravitational galaxy simulation implemented in Python and Pygame.

## Overview

The simulation models the motion of stars under gravitational forces.
The galaxy consists of a stellar disk and a central bulge, with stars
interacting through Newtonian gravity.

## Versions

### 1. Potential-based simulation

The initial version models the galaxy using a predefined gravitational
potential. Initial velocities are assigned to stars to produce an
approximately stable rotating disk.

### 2. Barnes-Hut algorithm

The final version uses a Barnes-Hut quadtree to approximate long-range
gravitational interactions.

Instead of computing the force between every pair of stars, distant
groups of stars are approximated by their combined mass and center of
mass, reducing the computational complexity from O(N²) to approximately
O(N log N).

## Features

- 2D gravitational N-body simulation
- Stellar disk and central bulge
- Initial orbital velocities
- Barnes-Hut quadtree
- Pygame visualization
- Configurable number of particles and simulation parameters

## Technologies

- Python
- NumPy
- Pygame
