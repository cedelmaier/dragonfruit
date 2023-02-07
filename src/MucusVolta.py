# General submission script for mucus simulations
#import hoomd
#import hoomd.md as md
#import gsd.hoomd

import argparse
import datetime
import itertools
import os
import sys
import yaml

import numpy as np

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), 'mucus'))
from mucus_seed import MucusSeed

# Create a rudimentary parser to get the number of steps, the write frequency,
# and the box length
def parse_args():
    parser = argparse.ArgumentParser(prog='MucusVolta.py')

    # General options
    parser.add_argument('--yaml', type=str, default='mucus.default.yaml',
            help='YAML configuration file')

    parser.add_argument('-d', '--workdir', action = 'store_true',
            help = 'Working directory')

    # Add verbosity control
    parser.add_argument('-v', '--verbose', action="store_true",
                        help = 'Verbose output')

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

    # Configurator is just a mucus seed
    configurator = MucusSeed(opts.workdir, opts)

    # Bail early if there is the LAMMPs simulation, just to make coding this easier for now
    # XXX There is probably a better way of doing this, but to make it use the same code, meh
    if configurator.engine == "LAMMPS":
        # Premature exit to not execute everything else
        configurator.Configure()
        configurator.PrintInformation()
        if opts.verbose: print(f"MucusVolta LAMMPS return")
        sys.exit(0)

    ###############################################################################
    # Set up the system
    ###############################################################################
    # Make sure that hoomd and gsd are imported
    import hoomd
    import hoomd.md as md
    import gsd.hoomd
    # Create the hoomd context on the CPU or GPU if appropriate
    if configurator.compute_mode == 'cpu':
        cpu = hoomd.device.CPU(notice_level = 2)
        sim = hoomd.Simulation(cpu, seed = configurator.nseed)
        snap = hoomd.Snapshot(cpu.communicator)
    elif configurator.compute_mode == 'gpu':
        gpu = hoomd.device.GPU(notice_level = 2)
        sim = hoomd.Simulation(gpu, seed = configurator.nseed)
        snap = hoomd.Snapshot(gpu.communicator)

    # Check how we are reading in information
    if configurator.init_type == 'create_equilibration':
        fudge_size = 0.2 # Prevent particles on the boundaries because reasons
    elif configurator.init_type == 'production':
        ftraj = gsd.hoomd.open(configurator.init_filename, 'rb')
        gsd_snap = ftraj[-1]
        # Convert to hoomd snapshot
        if configurator.compute_mode == 'cpu':
            snap = snap.from_gsd_snapshot(gsd_snap, cpu.communicator)
        elif configurator.compute_mode == 'gpu':
            snap = snap.from_gsd_snapshot(gsd_snap, gpu.communicator)
    else:
        print(f"Configurator init {configurator.init_type} not yet available, exiting!")
        sys.exit(1)

    # Run through the configuration and stuff
    configurator.Configure(snap)
    configurator.PrintInformation(snap)

    ###############################################################################
    # Create the system
    ###############################################################################
    sim.create_state_from_snapshot(snap)

    # Create the pair potentials for the system
    cell = hoomd.md.nlist.Cell(buffer = 0.4, exclusions = ['bond'])

    # If we are initializing, different stuff to do than if we are running 'production'
    if configurator.init_type == 'create_equilibration':
        # We have to ramp the potential by hand, which is annoying, so this is a completely different
        # sequence to do this
        # XXX Check if ramping or the soft potential is implemented in hoomd, as right now this sucks to do
        #v_prefactor = hoomd.variant.Ramp(0.0, 100.0, 0, configurator.nsteps_equilibrate)

        # Set up a default gauss potential so we can use it
        gauss = md.pair.Gauss(nlist=cell, default_r_cut=3.5, default_r_on=0.0)
        gauss.params[('E', 'E')] = {'epsilon': 0.0, 'sigma': 1.0}
        gauss.params[('E', 'C')] = {'epsilon': 0.0, 'sigma': 1.0}
        gauss.params[('E', 'H')] = {'epsilon': 0.0, 'sigma': 1.0}
        gauss.params[('E', 'P')] = {'epsilon': 0.0, 'sigma': 1.0}
        gauss.params[('C', 'C')] = {'epsilon': 0.0, 'sigma': 1.0}
        gauss.params[('C', 'H')] = {'epsilon': 0.0, 'sigma': 1.0}
        gauss.params[('C', 'P')] = {'epsilon': 0.0, 'sigma': 1.0}
        gauss.params[('H', 'H')] = {'epsilon': 0.0, 'sigma': 1.0}
        gauss.params[('H', 'P')] = {'epsilon': 0.0, 'sigma': 1.0}
        gauss.params[('P', 'P')] = {'epsilon': 0.0, 'sigma': 1.0}


    else:
        # WCA potential for the excluded volume interactions
        wca = md.pair.LJ(nlist=cell)
        wca.mode='shift'

        # LJ interaction for the attractive EM stuff
        lje = md.pair.LJ(nlist=cell)

        # Born mayer huggins potential for hydrophobic interactions
        bmh = md.pair.BornMayerHuggins(nlist=cell)

        ###############################
        # Excluded volume interactions
        ###############################
        # Set r_cut = 0 to disable interactions
        # Internal switches for disabling the lennard_jones and BMH interactions
        print(f"Detected WCA interaction for all <--> all")
        wca.params[('E', 'E')] = {'epsilon': 1.0, 'sigma': 1.0}
        wca.params[('E', 'C')] = {'epsilon': 1.0, 'sigma': 1.0}
        wca.params[('E', 'H')] = {'epsilon': 1.0, 'sigma': 1.0}
        wca.params[('E', 'P')] = {'epsilon': 1.0, 'sigma': 1.0} # Disabled?
        wca.params[('C', 'C')] = {'epsilon': 1.0, 'sigma': 1.0}
        wca.params[('C', 'H')] = {'epsilon': 1.0, 'sigma': 1.0}
        wca.params[('C', 'P')] = {'epsilon': 1.0, 'sigma': 1.0}
        wca.params[('H', 'H')] = {'epsilon': 1.0, 'sigma': 1.0} # Disabled?
        wca.params[('H', 'P')] = {'epsilon': 1.0, 'sigma': 1.0}
        wca.params[('P', 'P')] = {'epsilon': 1.0, 'sigma': 1.0}
        wca.r_cut[('E', 'E')] = np.power(2.0, 1/6)*1.0
        wca.r_cut[('E', 'C')] = np.power(2.0, 1/6)*1.0
        wca.r_cut[('E', 'H')] = np.power(2.0, 1/6)*1.0
        # Disable WCA if Lennard Jones present
        if configurator.lennard_jones <= 0.0:
            wca.r_cut[('E', 'P')] = np.power(2.0, 1/6*1.0)
        else:
            wca.r_cut[('E', 'P')] = 0.0
        wca.r_cut[('C', 'C')] = np.power(2.0, 1/6)*1.0
        wca.r_cut[('C', 'H')] = np.power(2.0, 1/6)*1.0
        wca.r_cut[('C', 'P')] = np.power(2.0, 1/6)*1.0
        # Disable WCA if BMH present
        if configurator.bmh >= 0.0:
            wca.r_cut[('H', 'H')] = np.power(2.0, 1/6)*1.0
        else:
            wca.r_cut[('H', 'H')] = 0.0
        wca.r_cut[('H', 'P')] = np.power(2.0, 1/6)*1.0
        wca.r_cut[('P', 'P')] = np.power(2.0, 1/6)*1.0

        ###############################
        # Attractive Lennard Jones type interactions
        ###############################
        # Only inlucde lennard jones if it is enabled
        if configurator.lennard_jones > 0.0:
            print(f"Detected Lennard-Jones attractive interaction for E <--> P")
            lje.params[('E', 'E')] = {'epsilon': configurator.lennard_jones, 'sigma': 1.0} # Disabled
            lje.params[('E', 'C')] = {'epsilon': configurator.lennard_jones, 'sigma': 1.0} # Disabled
            lje.params[('E', 'H')] = {'epsilon': configurator.lennard_jones, 'sigma': 1.0} # Disabled
            lje.params[('E', 'P')] = {'epsilon': configurator.lennard_jones, 'sigma': 1.0} 
            lje.params[('C', 'C')] = {'epsilon': configurator.lennard_jones, 'sigma': 1.0} # Disabled
            lje.params[('C', 'H')] = {'epsilon': configurator.lennard_jones, 'sigma': 1.0} # Disabled
            lje.params[('C', 'P')] = {'epsilon': configurator.lennard_jones, 'sigma': 1.0} # Disabled
            lje.params[('H', 'H')] = {'epsilon': configurator.lennard_jones, 'sigma': 1.0} # Disabled
            lje.params[('H', 'P')] = {'epsilon': configurator.lennard_jones, 'sigma': 1.0} # Disabled
            lje.params[('P', 'P')] = {'epsilon': configurator.lennard_jones, 'sigma': 1.0} # Disabled
            lje.r_cut[('E', 'E')] = 0.0
            lje.r_cut[('E', 'C')] = 0.0
            lje.r_cut[('E', 'H')] = 0.0
            lje.r_cut[('E', 'P')] = 2.5
            lje.r_cut[('C', 'C')] = 0.0
            lje.r_cut[('C', 'H')] = 0.0
            lje.r_cut[('C', 'P')] = 0.0
            lje.r_cut[('H', 'H')] = 0.0
            lje.r_cut[('H', 'P')] = 0.0
            lje.r_cut[('P', 'P')] = 0.0

        ###############################
        # Attractive Born-Mayer-Huggins interactions
        ###############################
        # Only inlucde if used (attractive)
        if configurator.bmh < 0.0:
            print(f"Detected Born-Mayer-Huggins attractive interaction for H <--> H")
            bmh.params[('E', 'E')] = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            bmh.params[('E', 'C')] = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            bmh.params[('E', 'H')] = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            bmh.params[('E', 'P')] = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            bmh.params[('C', 'C')] = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            bmh.params[('C', 'H')] = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            bmh.params[('C', 'P')] = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            bmh.params[('H', 'H')] = {'A': configurator.bmh, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            bmh.params[('H', 'P')] = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            bmh.params[('P', 'P')] = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            bmh.r_cut[('E', 'E')] = 0.0
            bmh.r_cut[('E', 'C')] = 0.0
            bmh.r_cut[('E', 'H')] = 0.0
            bmh.r_cut[('E', 'P')] = 0.0
            bmh.r_cut[('C', 'C')] = 0.0
            bmh.r_cut[('C', 'H')] = 0.0
            bmh.r_cut[('C', 'P')] = 0.0
            bmh.r_cut[('H', 'H')] = 2.5
            bmh.r_cut[('H', 'P')] = 0.0
            bmh.r_cut[('P', 'P')] = 0.0
        

    ###############################
    # Bonded and angle interactions
    ###############################
    # Assign bonded interaction strengths
    fenewca = md.bond.FENEWCA()
    fenewca.params['mucusbond'] = dict(k = 30.0, r0 = 1.5, epsilon = 1.0, sigma = 1.0, delta = 0.0)

    # Only add angle if there is an angular potential, otherwise wasting space
    if configurator.mucusbend > 0.0:
        angleharmonic = md.angle.Harmonic()
        angleharmonic.params['mucusbend'] = dict(k = 0.0, t0 = np.pi)

    # Set up the integration for the system
    integrator = md.Integrator(dt = configurator.deltatau)
    langevin = md.methods.Langevin(hoomd.filter.All(),
                                   kT = configurator.kT)
    # Gamma for every particle is mass/t_damp
    langevin.gamma['E'] = 1.0/configurator.t_damp
    langevin.gamma['C'] = 1.0/configurator.t_damp
    langevin.gamma['H'] = 1.0/configurator.t_damp
    langevin.gamma['P'] = 1.0/configurator.t_damp

    integrator.methods.append(langevin)

    # Append the forces that we need
    if configurator.init_type == 'create_equilibration':
        integrator.forces.append(gauss)
    elif configurator.init_type == 'production':
        integrator.forces.append(wca)
        if configurator.lennard_jones > 0.0:
            integrator.forces.append(lje)
        if configurator.bmh < 0.0:
            integrator.forces.append(bmh)
    integrator.forces.append(fenewca)
    if configurator.mucusbend > 0.0:
        integrator.forces.append(angleharmonic)

    # Add to the simulation
    sim.operations.integrator = integrator

    # Run a simulation to get variables tow ork out, then thermalize the system
    #pre_snapshot = sim.state.get_snapshot()
    #print(pre_snapshot)
    #print(pre_snapshot.particles)
    #print(pre_snapshot.particles.position[144])
    #print(pre_snapshot.particles.velocity[144])
    #print(pre_snapshot.particles.acceleration[144])
    #for ibond in range(len(pre_snapshot.bonds.typeid)):
    #    print(f"  Bond[{ibond}] type {pre_snapshot.bonds.typeid[ibond]} group {pre_snapshot.bonds.group[ibond]}")
    #sys.exit(1)

    if configurator.init_type == 'create_equilibration':
        sim.run(0)
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
    if configurator.init_type == 'create_equilibration':
        logger.add(gauss, quantities=['energies', 'forces'])
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

    # Set up writing out GSD trajectories
    gsd_writer = None
    if configurator.init_type == 'create_equilibration':
        gsd_writer = hoomd.write.GSD(filename = configurator.trajectory_file,
                                     trigger = hoomd.trigger.Periodic(configurator.nwrite_equilibrate),
                                     mode = 'wb',
                                     filter = hoomd.filter.All(),
                                     log = logger)
    else:
        gsd_writer = hoomd.write.GSD(filename = configurator.trajectory_file,
                                     trigger = hoomd.trigger.Periodic(configurator.nwrite),
                                     mode = 'wb',
                                     filter = hoomd.filter.All(),
                                     log = logger)
    sim.operations.writers.append(gsd_writer)

    print(f"--------")
    if configurator.init_type == 'create_equilibration':
        print(f"Equilibrating...")

        # Divide this up into 1000 different regimes for nsteps_equilibrate
        nblocks = 1000
        nblock_size = np.int32(configurator.nsteps_equilibrate / nblocks)
        for iblock in range(1000):
            print(f"  Block {iblock}")
            gauss.params[('E', 'E')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
            gauss.params[('E', 'C')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
            gauss.params[('E', 'H')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
            gauss.params[('E', 'P')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
            gauss.params[('C', 'C')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
            gauss.params[('C', 'H')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
            gauss.params[('C', 'P')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
            gauss.params[('H', 'H')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
            gauss.params[('H', 'P')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
            gauss.params[('P', 'P')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
            sim.run(nblock_size)

    elif configurator.init_type == 'production':
        sim.run(configurator.nsteps)

