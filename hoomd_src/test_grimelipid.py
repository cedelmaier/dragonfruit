import itertools
import math
import sys

import gsd.hoomd
import hoomd
import numpy

import os

# Compute mode for this test
compute_mode = 'cpu'
lbox = 100

# Switch on compute mode to check what is going on
cpu = None
sim = None
snap = None
if compute_mode == 'cpu':
    cpu = hoomd.device.CPU(notice_level = 2)
    sim = hoomd.Simulation(device = cpu, seed = 1)
    snap = hoomd.Snapshot(cpu.communicator)

snap.particles.N = 2
snap.particles.position[0] = [0, 0, 0]
snap.particles.position[1] = [0.0, 0.0, 0.75]
snap.particles.typeid[0] = 0
snap.particles.typeid[1] = 0

snap.configuration.box = [lbox, lbox, lbox, 0, 0, 0]

snap.particles.types = ['A']

# Create the snapshot state
sim.create_state_from_snapshot(snap)

# Create a neighbor list, potential, all that jazz
nl = hoomd.md.nlist.Cell(buffer = 0.4)
glf = hoomd.md.pair.GrimeLipid(nlist = nl)
glf.params[('A', 'A')] = {'A': 10.0, 'B': 2.0, 'r0': 0.6, 'rc': 1.2}
glf.r_cut[('A', 'A')] = 2.5

# Creat the integrator, etc
integrator = hoomd.md.Integrator(dt = 0.005)
langevin = hoomd.md.methods.Langevin(filter = hoomd.filter.All(), kT = 1.0)
integrator.methods.append(langevin)
integrator.forces.append(glf)

sim.operations.integrator = integrator
sim.run(0)

# Get the information back out of the simulation
snap_res = sim.state.get_snapshot()
sim_energies = sim.operations.integrator.forces[0].energies
sim_forces = sim.operations.integrator.forces[0].forces

if snap_res.communicator.rank == 0:
    print(f"Energies: {sim_energies}")
    print(f"Forces: {sim_forces}")

