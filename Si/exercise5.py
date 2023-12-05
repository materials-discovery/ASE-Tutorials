from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.optimize.bfgs import BFGS
from ase.constraints import ExpCellFilter
# Create a Silicon unit cell
a = 5.43  # lattice constant for Silicon in Angstrom
si = bulk('Si', 'diamond', a=a)

# Set up Quantum ESPRESSO calculator

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
                tstress=True, tprnfor=True, kpts=('KPOINT', 'KPOINT', 'KPOINT'), input_data=input_data)

# Attach the calculator to the Silicon atoms
si.calc = calc

# Perform a geometry optimization

stress_constraints = ExpCellFilter(si)

opt = BFGS(stress_constraints, trajectory='Opt.traj', logfile='Opt.log')
opt.run(fmax=0.01)
