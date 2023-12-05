from ase.calculators.espresso import Espresso
from ase.io import read
from ase.dft.kpoints import BandPath
from ase.spectrum.band_structure import BandStructure
from ase.dft import DOS 
import numpy as np 
import matplotlib.pyplot as plt 

si = read('Opt.traj@-1')

pseudopotentials = {'Si': 'Si.vbc.UPF'}

input_data = {
    'control': {
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
        'prefix': 'si_scf',
        'outdir': './',
    },
    'system': {
        'ecutwfc': 50,
        'ecutrho': 200,
        'occupations': 'smearing',
        'smearing': 'mv',
        'degauss': 0.01,
        'input_dft': 'PBE',
    },
    'electrons': {
        'conv_thr': 1.0e-8,
    },
}

# Set up the Quantum ESPRESSO calculator

calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=(6, 6, 6), input_data=input_data)

# Attach the calculator to the Silicon atoms
si.calc = calc 

si.get_potential_energy()


# Get the Fermi level
fermi_level = calc.get_fermi_level()

# Get the lattice and band path
lattice = si.get_cell()
bandpath = lattice.bandpath(npoints=100)

# Set up the calculator for band structure calculation
input_data['control'].update({'calculation':'bands','restart_mode':'restart','verbosity':'high'})

calc.set(kpts=bandpath, input_data=input_data)
calc.calculate(si)

# Get the band structure and modify the energies
bs = calc.band_structure()
new_band = BandStructure(energies=bs.energies - fermi_level * np.ones_like(bs.energies), path=bs.path)

new_band.plot()
# Plot and save the band structure

plt.savefig("band.png", dpi=300)

