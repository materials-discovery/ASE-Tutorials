# Exercise 1: Demonstrates molecular dynamics with constant energy.

from ase import units
from ase.build import bulk, make_supercell
from ase.calculators.emt import EMT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.io.trajectory import Trajectory
from ase.optimize.bfgs import BFGS
from ase.constraints import StrainFilter
import numpy as np
from ase.visualize import view


# Set up a crystal
element_symbol = 'Cu'
atoms = bulk(element_symbol, 'fcc', a=3.597)

# Describe the interatomic interactions with the Effective Medium Theory
atoms.calc = EMT()

# Before we run the molecular dynamics, let's first optimize the unit cell
constraints = StrainFilter(atoms)

trajectory_filename = f'{element_symbol}_opt.traj'
logfile_filename = f'{element_symbol}_opt.log'

opt = BFGS(constraints, trajectory=trajectory_filename, logfile=logfile_filename)
opt.run(fmax=0.01)

# Let's view the atoms
view(atoms, viewer='ngl')

# After optimization, access the optimized lattice constant 'a'
optimized_a = atoms.get_cell_lengths_and_angles()[0] * (2 ** 0.5)  # For a cubic cell, the first value represents 'a'

print(f"The optimized lattice constant 'a' is: {optimized_a} Å")


# Calculate the volume of the unit cell for FCC
volume_unit_cell = optimized_a ** 3

# For FCC, there are 4 atoms per unit cell
num_atoms_unit_cell = 4

# FCC packing density calculation
radius = optimized_a / (2 * (2 ** 0.5))  # Radius of the atoms in an FCC lattice

# Calculate the volume of one sphere in the unit cell
volume_one_sphere = (4 / 3) * np.pi * (radius ** 3)

# Calculate the packing density
packing_density = ((num_atoms_unit_cell * volume_one_sphere) / volume_unit_cell) * 100

print(f"The packing efficiency of the optimized FCC unit cell is: {packing_density:.4f}%")


# Now, let's create a supercell for molecular dynamics
supercell_size = [[5, 0, 0], [0, 5, 0], [0, 0, 5]]  # Define the size of the supercell (5x5x5)
supercell = make_supercell(atoms, supercell_size)

# Let's view the supercell
view(supercell, viewer='ngl')

supercell.calc = EMT()

# Set the momenta corresponding to T=300K
MaxwellBoltzmannDistribution(supercell, temperature_K=300)

# Run MD with constant energy using the VelocityVerlet algorithm.
dyn = VelocityVerlet(supercell, 5 * units.fs, logfile='Cu_md.log')  # 5 fs time step.

def printenergy(a=supercell):
    """Function to print the potential, kinetic, and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

# Run the dynamics
dyn.attach(printenergy, interval=10)

trajectory_filename = f'{element_symbol}_md.traj'
logfile_filename = f'{element_symbol}_md.log'

traj = Trajectory(trajectory_filename, 'w', supercell)
dyn.attach(traj.write, interval=10)

printenergy()
dyn.run(200)


