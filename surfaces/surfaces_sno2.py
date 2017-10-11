#!/usr/bin/env python

from __future__ import print_function

import os, sys, json

import numpy as np

from ase.build import bulk, surface
from ase.units import Rydberg, Bohr
from ase.io import read
from ase.visualize import view
from ase.spacegroup import crystal 
from ase.calculators.espresso import Espresso

infile = sys.argv[1]
print(infile)
with open(infile) as handle:
    system = json.loads(handle.read())

face = system['face']
layers = system['layers']
vacuum = system['vacuum']
kpts = system['kpts']
ecut = system['ecut']
mode = system['mode']

cwd = os.getcwd()

def cassiterite(show=False):
    a = 4.7382
    c = 3.1871 
    sno2 = crystal(['Sn','O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                   spacegroup='P 4 2/mnm', cellpar=[a,a,c,90,90,90], 
                   pbc=True)
    return sno2


def create_surface(atoms,face=(1,1,0),layers=3,vacuum=10.813,kpts=([6,3,1])):
    mySurface = surface( atoms, face, layers)
    mySurface.center(vacuum=vacuum, axis=2)
    kpts = np.asarray(kpts)
    return mySurface


sno2         = cassiterite()
sno2_surface = create_surface(sno2,
                              face=face,
                              layers=layers,
                              vacuum=vacuum,
                              kpts=kpts)

# Put together QE input dict
input_dict = {
    'control': {
        'calculation': 'scf',
        'etot_conv_thr': 1e-6,
        'nstep': 100,
        'outdir': 'sno2_test_face_{0}{1}{2}'.format(face[0], face[1], face[2]),
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


sno2_surface.set_calculator(calc)

if mode == 'view':
    view(sno2_surface)
elif mode == 'calc':
    calc.calculate(sno2_surface)
    print('SnO2 PE:', sno2_surface.get_potential_energy())

