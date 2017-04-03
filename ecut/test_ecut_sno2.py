#!/usr/bin/env python

from __future__ import print_function

import os, sys, json

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

minecut = system['min_ecut']
maxecut = system['max_ecut']
intecut = system['interval']
kpts = system['kpts']
mode = system['mode']

cwd = os.getcwd()

def cassiterite(show=False):
    a = 4.7382
    c = 3.1871
    sno2 = crystal(['Sn','O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                   spacegroup='P 4 2/mnm', cellpar=[a,a,c,90,90,90], 
                   pbc=True)
    return sno2


if mode == "view":
    sno2 = cassiterite()
    view(sno2)

elif mode == "calc":
    ecut_dict = dict()
    for ecut in range(minecut,maxecut+intecut,intecut):
        print('Working on...',ecut,' -> ', end=' ')
        sno2 = cassiterite()
        calc = Espresso(pw=ecut * Rydberg, calculation='scf', kpts=kpts,
                        psppath=cwd+"/../pseudo",
                        convergence={'energy': 1e-6,
                                     'maxsteps': 100, 'diag': 'cg'},
                        outdir='sno2_test',
                        site='local'
                        )
        sno2.set_calculator(calc)
        calc.calculate(sno2)
        ecut_dict[ecut] = sno2.get_potential_energy()
        print('ECUT:', ecut, 'SnO2 PE:', sno2.get_potential_energy())
        
    emin = min(ecut_dict.itervalues())
    for key in ecut_dict:
        ecut_dict[key] -= emin

    import matplotlib.pyplot as plt
    x, y = zip( *sorted(ecut_dict.items()) )
    plt.plot(x,y,'-o')
    plt.xlabel(r'$E_{cut} (Ry) $')
    plt.ylabel(r'$\Delta E (eV) $')
    plt.show()


