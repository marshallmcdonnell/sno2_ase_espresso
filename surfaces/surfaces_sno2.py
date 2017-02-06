#!/usr/bin/env python

from __future__ import print_function

import sys, json

import numpy as np

from ase.build import bulk, surface
from ase.units import Rydberg, Bohr
from ase.io import read
from ase.visualize import view
from ase.spacegroup import crystal 
from espresso import Espresso

infile = sys.argv[1]
print(infile)
with open(infile) as handle:
    system = json.loads(handle.read())

face = system['face']
layers = system['layers']
vacuum = system['vacuum']
kpts = system['kpts']
ecut = system['ecut']


def cassiterite(show=False):
    a = 4.7382
    c = 3.1871 
    sno2 = crystal(['Sn','O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                   spacegroup='P 4 2/mnm', cellpar=[a,a,c,90,90,90], 
                   pbc=True)
    return sno2


def create_surface(sno2,face=(1,1,0),layers=3,vacuum=10.813,kpts=([6,3,1])):
    sno2_surface = surface( sno2, face, layers)
    sno2_surface.center(vacuum=vacuum, axis=2)
    kpts = np.asarray(kpts)
    return sno2_surface


sno2         = cassiterite()
sno2_surface = create_surface(sno2,
                              face=face,
                              layers=layers,
                              vacuum=vacuum,
                              kpts=kpts)

calc = Espresso(pw=ecut * Rydberg, calculation='scf', kpts=kpts,
                convergence={'energy': 1e-6,
                             'maxsteps': 100, 'diag': 'cg'},
                psppath="/home/ntm/projects/josh_kim/sno2/pseudo",
                outdir='sno2_test'
                )

#view(sno2_surface)
sno2_surface.set_calculator(calc)

calc.calculate(sno2_surface)
print('SnO2 PE:', sno2_surface.get_potential_energy())

