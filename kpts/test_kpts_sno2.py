#!/usr/bin/env python

from __future__ import print_function

import numpy as np

from ase.build import bulk, surface
from ase.units import Rydberg, Bohr
from ase.io import read
from ase.visualize import view
from ase.spacegroup import crystal 
from espresso import Espresso

def cassiterite(show=False):
    a = 4.7382
    c = 3.1871
    sno2 = crystal(['Sn','O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                   spacegroup='P 4 2/mnm', cellpar=[a,a,c,90,90,90], 
                   pbc=True)
    return sno2


ecut=35.0
kpts_dict = dict()
for k in range(1,7):
    print('Working on...',k,' -> ', end=' ')
    kpts = (k,k,k)
    sno2 = cassiterite()
    calc = Espresso(pw=ecut * Rydberg, calculation='scf', kpts=kpts,
                    psppath="/home/ntm/projects/josh_kim/sno2/pseudo",
                    convergence={'energy': 1e-6,
                                 'maxsteps': 100, 'diag': 'cg'},
                    outdir='sno2_test'
                    )
    sno2.set_calculator(calc)
    calc.calculate(sno2)
    kpts_dict[k] = sno2.get_potential_energy()
    print('K_POINTS:', k, 'SnO2 PE:', sno2.get_potential_energy())
    

import matplotlib.pyplot as plt
x, y = zip( *sorted(kpts_dict.items()) )
plt.plot(x,y,'-o')
plt.show()


