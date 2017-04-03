#!/usr/bin/env python

from __future__ import print_function

import os, sys,json

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

min_kpts = system['min_kpts']
max_kpts = system['max_kpts']
interval = system['interval']
ecut = system['ecut']
mode = system['mode']

cwd = os.getcwd()

kpts_dict = dict()
def cassiterite(show=False):
    a = 4.7382
    c = 3.1871
    sno2 = crystal(['Sn','O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                   spacegroup='P 4 2/mnm', cellpar=[a,a,c,90,90,90], 
                   pbc=True)
    return sno2

if mode == 'view':
    sno2 = cassiterite()
    view(sno2)

if mode == 'calc':
    for k in range(min_kpts,max_kpts+interval,interval):
        print('Working on...',k,' -> ', end=' ')
        kpts = (k,k,k)
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
        kpts_dict[k] = sno2.get_potential_energy()
        print('K_POINTS:', k, 'SnO2 PE:', sno2.get_potential_energy())

    emin = min(kpts_dict.itervalues())
    for key in kpts_dict:
        kpts_dict[key] -= emin

        

    import matplotlib.pyplot as plt
    x, y = zip( *sorted(kpts_dict.items()) )
    plt.plot(x,y,'-o')
    plt.xlabel(r'K-point mesh dimensions (cubic)')
    plt.ylabel(r'$\Delta E (eV) $')
    plt.show()

    plt.show()


