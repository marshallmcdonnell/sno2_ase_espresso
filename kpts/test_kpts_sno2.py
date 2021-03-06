#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import json

from ase.visualize import view
from ase.spacegroup import crystal
from ase.calculators.espresso import Espresso

# Read input JSON file
infile = sys.argv[1]
print(infile)
with open(infile) as handle:
    system = json.loads(handle.read())

min_kpts = system['min_kpts']
max_kpts = system['max_kpts']
interval = system['interval']
ecut = system['ecut']
mode = system['mode']


# Structure for SnO2
def cassiterite(show=False):
    a = 4.7382
    c = 3.1871
    sno2 = crystal(['Sn', 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                   spacegroup='P 4 2/mnm', cellpar=[a, a, c, 90, 90, 90],
                   pbc=True)
    return sno2

# Just view the structure generated
if mode == "view":
    sno2 = cassiterite()
    view(sno2)

# Actually run the calculation
elif mode == "calc":

    # Dict to save output from each calculation
    kpts_dict = dict()

    # Save current working directory 
    cwd = os.getcwd()

    # Loop over k-point mesh to check for convergence
    for k in range(min_kpts,max_kpts+interval,interval):
        print('Working on...', k, ' -> ', end=' ')

        # k-point mesh
        kpts = (k,k,k)

        # Make the structure
        sno2 = cassiterite()

        # Put together QE input dict
        input_dict = {
            'control': {
                'calculation': 'scf',
                'etot_conv_thr': 1e-6,
                'nstep': 100,
                'outdir': 'sno2_kpts_{0}x{1}x{2}'.format(k,k,k),
            },
            'system': {
                'ecutwfc': ecut,
            },
            'electrons': {
                'diagonalization': 'cg',
            },
        }

        # Put together pseudopotential dict
        psp_dict = {'Sn': 'Sn.UPF',
                    'O': 'O.UPF',
                    }

        calc = Espresso(input_data=input_dict,
                        kpts=kpts,
                        pseudo_dir=cwd + "/../pseudo",
                        pseudopotentials=psp_dict,
                        )

        # Attach calculator to structure
        sno2.set_calculator(calc)

        # Calculate and save output
        calc.calculate(sno2)
        kpts_dict[k] = sno2.get_potential_energy()
        print('K_POINTS:', k, 'SnO2 PE:', sno2.get_potential_energy())

    # Get the energy minimum and adjust to the '0' value to get delta_E
    emin = min(kpts_dict.itervalues())
    for key in kpts_dict:
        kpts_dict[key] -= emin

    # Plot result
    import matplotlib.pyplot as plt
    x, y = zip( *sorted(kpts_dict.items()) )
    plt.plot(x,y,'-o')
    plt.xlabel(r'K-point mesh dimensions (cubic)')
    plt.ylabel(r'$\Delta E (eV) $')
    plt.show()
