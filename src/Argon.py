# Put a licencse here or something

# This is based off of LipidVolta.py, but should be for an easier system
import hoomd
import hoomd.md as md
import gsd.hoomd

import argparse
import os
import sys
import yaml

import numpy as np

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), 'argon'))
from argon_seed import ArgonSeed

# Create a rudimentary parser to get the number of steps, the write frequency,
# and the box length
def parse_args():
    parser = argparse.ArgumentParser(prog='Argon.py')

    # General options
    parser.add_argument('--yaml', type=str, default='lipid.default.yaml',
            help='YAML configuration file')

    parser.add_argument('-d', '--workdir', action = 'store_true',
            help = 'Working directory')


    opts = parser.parse_args()

    return opts

# Create a status line maker for our output
class Status():

    def __init__(self, sim):
        self.sim = sim

    @property
    def seconds_remaining(self):
        try:
            return (self.sim.final_timestep - self.sim.timestep) / self.sim.tps
        except ZeroDivisionError:
            return 0

    @property
    def etr(self):
        return str(datetime.timedelta(seconds=self.seconds_remaining))

###############################################################################
# Main program start
###############################################################################
if __name__ == "__main__":
    # Parse those args
    opts = parse_args()

    cwd = os.getcwd()
    if not opts.workdir:
        opts.workdir = os.path.abspath(cwd)
    elif not os.path.exists(opts.workdir):
        raise IOError("Working directory {} does not exist, exiting.".format(
            opts.workdir))
    else:
        opts.workdir = os.path.abspath(opts.workdir)

    # Configurator is just a septin seed, since we moved the functionality into that class
    configurator = ArgonSeed(opts.workdir, opts)

    ###############################################################################
    # Set up the system
    ###############################################################################
    # Create the hoomd context on the GPU
    if configurator.compute_mode == 'cpu':
        cpu = hoomd.device.CPU(notice_level = 2)
        sim = hoomd.Simulation(cpu, seed = configurator.nseed)
        snap = hoomd.Snapshot(cpu.communicator)
    elif configurator.compute_mode == 'gpu':
        gpu = hoomd.device.GPU(notice_level = 2)
        sim = hoomd.Simulation(gpu, seed = configurator.nseed)
        snap = hoomd.Snapshot(gpu.communicator)

    # Initialize the argon simulation
    configurator.Configure(snap)

    # Print the information for the simulation
    configurator.PrintInformation(snap)

    ###############################################################################
    # Create the system
    ###############################################################################
    sim.create_state_from_snapshot(snap)

    # Create the pair potential
    cell = hoomd.md.nlist.Cell(buffer=0.4)
    # XXX: Figure out how to create a simple LJ potential here!
    # Somethig like lj = ...?

    # Set up the integrator for the system
    # XXX: Figure out how to set up the integrator for the system in the NVE ensemble for liquid argon

    # Run a timestep to get variables and then thermalize
    sim.state.thermalize_particle_momenta(hoomd.filter.All(), configurator.kT)

    ###############################################################################
    # Print information for the main program
    ###############################################################################
    # Keep track of the thermodynamic information
    thermodynamic_properties = md.compute.ThermodynamicQuantities(
            filter = hoomd.filter.All())
    sim.operations.computes.append(thermodynamic_properties)
    
    # Set up the logging of important quantities
    logger = hoomd.logging.Logger()
    logger.add(lj, quantities=['energies', 'forces'])
    logger.add(sim, quantities=['timestep', 'walltime', 'tps'])
    logger.add(thermodynamic_properties)
    
    # Display some quantities to a table while running
    output_logger = hoomd.logging.Logger(categories=['scalar', 'string'])
    status = Status(sim)
    output_logger.add(sim, quantities=['timestep', 'tps'])
    output_logger[('Status', 'etr')] = (status, 'etr', 'string')
    output_logger.add(thermodynamic_properties, quantities=['kinetic_temperature', 'pressure'])
    table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(period=configurator.nwrite),
                              logger=output_logger)
    sim.operations.writers.append(table)

    # Set up writing out to a GSD file for trajectories
    gsd_writer = hoomd.write.GSD(filename = configurator.trajectory_file,
                                 trigger = hoomd.trigger.Periodic(configurator.nwrite),
                                 mode = 'wb',
                                 filter = hoomd.filter.All(),
                                 log = logger)
    sim.operations.writers.append(gsd_writer)
    
    print(f"--------")
    print(f"Running...")
    sim.run(configurator.nsteps)
    
