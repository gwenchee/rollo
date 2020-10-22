import openmc
import numpy as np
from numpy import sin, cos, tan, pi

# Constants
T_r1 = 2135e-5
T_r2 = 3135e-5
T_r3 = 3485e-5
T_r4 = 3835e-5
T_r5 = 4235e-5

uoc_9 = openmc.Material()
uoc_9.set_density('g/cc', 11)
uoc_9.add_nuclide('U235', 2.27325e-3)
uoc_9.add_nuclide('U238', 2.269476e-2)
uoc_9.add_nuclide('O16', 3.561871e-2)
uoc_9.add_nuclide('C0', 9.79714e-3)
uoc_9.temperature = 1110
uoc_9.volume = 4 / 3 * pi * (T_r1 ** 3) * 101 * 210 * 4 * 36

por_c = openmc.Material()
por_c.set_density('g/cc', 1)
por_c.add_nuclide('C0', 5.013980e-2)
por_c.temperature = 948

si_c = openmc.Material()
si_c.set_density('g/cc', 3.2)
si_c.add_nuclide('Si28', 4.431240e-2)
si_c.add_nuclide('Si29', 2.25887e-3)
si_c.add_nuclide('Si30', 1.48990e-3)
si_c.add_nuclide('C0', 4.806117e-2)
si_c.temperature = 948

graphite = openmc.Material()
graphite.set_density('g/cc', 1.8)
graphite.add_nuclide('C0', 9.025164e-2)
graphite.temperature = 948

lm_graphite = openmc.Material()
lm_graphite.set_density('g/cc', 1.8)
lm_graphite.add_nuclide('C0', 9.025164e-2)
lm_graphite.temperature = 948

flibe = openmc.Material()
flibe.set_density('g/cc', 1.95)
flibe.add_nuclide('Li6', 1.383014e-6)
flibe.add_nuclide('Li7', 2.37132e-2)
flibe.add_nuclide('Be9', 1.18573e-2)
flibe.add_nuclide('F19', 4.74291e-2)
flibe.temperature = 948

mats = openmc.Materials(
    (uoc_9,
     por_c,
     si_c,
     graphite,
     lm_graphite,
     flibe))

spheres = [openmc.Sphere(r=r)
           for r in [T_r1, T_r2, T_r3, T_r4, T_r5]]
triso_cells = [openmc.Cell(fill=uoc_9, region=-spheres[0]),
               openmc.Cell(fill=por_c, region=+spheres[0] & -spheres[1]),
               openmc.Cell(fill=graphite, region=+spheres[1] & -spheres[2]),
               openmc.Cell(fill=si_c, region=+spheres[2] & -spheres[3]),
               openmc.Cell(fill=graphite, region=+spheres[3] & -spheres[4])]
triso_univ = openmc.Universe(cells=triso_cells)

