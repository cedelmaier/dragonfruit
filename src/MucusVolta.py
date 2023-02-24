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

    # Set up a common adjusted radius if we have size differences
    r_sphere    = 0.5
    r_histone   = configurator.r_histone
    adj_diameter = 0.5 + configurator.r_histone

    # Figure out which (and how many) neighbor lists to start up
    print(f"Creating neighbor lists")
    # Create a variable number no matter what
    nlists = []
    for ilist in range(configurator.nlist_n):
        print(f"  Creating neighbor list {ilist}: {configurator.nlist_type[ilist]}, buffer: {configurator.nlist_buffer[ilist]}")
        nl = None
        if configurator.nlist_type[ilist] == 'cell':
            nl = hoomd.md.nlist.Cell(buffer = configurator.nlist_buffer[ilist], exclusions = ['bond'])
        elif configurator.nlist_type[ilist] == 'tree':
            nl = hoomd.md.nlist.Tree(buffer = configurator.nlist_buffer[ilist], exclusions = ['bond'])
        else:
            print(f"ERROR: Neighbor list type {configurator.nlist_type[ilist]} not currently supported, exiting!")
            sys.exit(1)
        nlists.append(nl)

    print(f"Setting interactions")
    # If we are initializing, different stuff to do than if we are running 'production'
    if configurator.init_type == 'create_equilibration':
        # We have to ramp the potential by hand, which is annoying, so this is a completely different
        # sequence to do this. Set up a default gaussian potential to use
        if configurator.equilibration_potential == 'gauss':
            print(f"  Equilibration: gaussian")
            gauss = md.pair.Gauss(nlist = nlists[0], default_r_cut = 3.5, default_r_on = 0.0)
            gauss.params[('muc_e', 'muc_e')] = {'epsilon': 0.0, 'sigma': 1.0}
            gauss.params[('muc_e', 'muc_c')] = {'epsilon': 0.0, 'sigma': 1.0}
            gauss.params[('muc_e', 'muc_h')] = {'epsilon': 0.0, 'sigma': 1.0}
            gauss.params[('muc_e', 'muc_p')] = {'epsilon': 0.0, 'sigma': 1.0*adj_diameter}
            gauss.params[('muc_c', 'muc_c')] = {'epsilon': 0.0, 'sigma': 1.0}
            gauss.params[('muc_c', 'muc_h')] = {'epsilon': 0.0, 'sigma': 1.0}
            gauss.params[('muc_c', 'muc_p')] = {'epsilon': 0.0, 'sigma': 1.0*adj_diameter}
            gauss.params[('muc_h', 'muc_h')] = {'epsilon': 0.0, 'sigma': 1.0}
            gauss.params[('muc_h', 'muc_p')] = {'epsilon': 0.0, 'sigma': 1.0*adj_diameter}
            gauss.params[('muc_p', 'muc_p')] = {'epsilon': 0.0, 'sigma': 2.0*configurator.r_histone}
            gauss.r_cut[('muc_e', 'muc_p')] = 3.5*adj_diameter
            gauss.r_cut[('muc_c', 'muc_p')] = 3.5*adj_diameter
            gauss.r_cut[('muc_h', 'muc_p')] = 3.5*adj_diameter
            gauss.r_cut[('muc_p', 'muc_p')] = 3.5*2.0*configurator.r_histone
        else:
            print(f"  Equilibration: GrimeLipid (custom)")
            # Use our own grime-lipid potential, as that has a soft intraction, and a defined cutoff...
            glf = md.pair.GrimeLipid(nlist = nlists[0], default_r_cut = 1.0)
            glf.params[('muc_e', 'muc_e')] = {'A': 0.0, 'B': 0.0, 'r0': 1.0, 'rc': 2.0}
            glf.params[('muc_e', 'muc_c')] = {'A': 0.0, 'B': 0.0, 'r0': 1.0, 'rc': 2.0}
            glf.params[('muc_e', 'muc_h')] = {'A': 0.0, 'B': 0.0, 'r0': 1.0, 'rc': 2.0}
            glf.params[('muc_e', 'muc_p')] = {'A': 0.0, 'B': 0.0, 'r0': adj_diameter, 'rc': 2.0*adj_diameter}
            glf.params[('muc_c', 'muc_c')] = {'A': 0.0, 'B': 0.0, 'r0': 1.0, 'rc': 2.0}
            glf.params[('muc_c', 'muc_h')] = {'A': 0.0, 'B': 0.0, 'r0': 1.0, 'rc': 2.0}
            glf.params[('muc_c', 'muc_p')] = {'A': 0.0, 'B': 0.0, 'r0': adj_diameter, 'rc': 2.0*adj_diameter}
            glf.params[('muc_h', 'muc_h')] = {'A': 0.0, 'B': 0.0, 'r0': 1.0, 'rc': 2.0}
            glf.params[('muc_h', 'muc_p')] = {'A': 0.0, 'B': 0.0, 'r0': adj_diameter, 'rc': 2.0*adj_diameter}
            glf.params[('muc_p', 'muc_p')] = {'A': 0.0, 'B': 0.0, 'r0': 2.0*r_histone, 'rc': 4.0*r_histone}

        # If we have a second neighbor list, we know that it goes to the muc_p group
        if configurator.nlist_n > 1:
            print(f"ERROR: Currently only 1 neighbor list for equilibration, exiting!")
            sys.exit(1)
            # Create a second interaction potential for just the muc_p <--> anything, and set the default r_cut to 0.0 to disable all interactions unless needed
            # Also create a default version of the parameter control to make hoomd happy
            gauss_large = md.pair.Gauss(nlist = nlists[1], default_r_cut = 0.0, default_r_on = 0.0)
            gauss_large.params[('muc_e', 'muc_e')] = {'epsilon': 0.0, 'sigma': 1.0}
            gauss_large.params[('muc_e', 'muc_c')] = {'epsilon': 0.0, 'sigma': 1.0}
            gauss_large.params[('muc_e', 'muc_h')] = {'epsilon': 0.0, 'sigma': 1.0}
            gauss_large.params[('muc_e', 'muc_p')] = {'epsilon': 0.0, 'sigma': 1.0*adj_diameter}
            gauss_large.params[('muc_c', 'muc_c')] = {'epsilon': 0.0, 'sigma': 1.0}
            gauss_large.params[('muc_c', 'muc_h')] = {'epsilon': 0.0, 'sigma': 1.0}
            gauss_large.params[('muc_c', 'muc_p')] = {'epsilon': 0.0, 'sigma': 1.0*adj_diameter}
            gauss_large.params[('muc_h', 'muc_h')] = {'epsilon': 0.0, 'sigma': 1.0}
            gauss_large.params[('muc_h', 'muc_p')] = {'epsilon': 0.0, 'sigma': 1.0*adj_diameter}
            gauss_large.params[('muc_p', 'muc_p')] = {'epsilon': 0.0, 'sigma': 2.0*configurator.r_histone}
            gauss_large.r_cut[('muc_e', 'muc_p')] = 3.5*adj_diameter
            gauss_large.r_cut[('muc_c', 'muc_p')] = 3.5*adj_diameter
            gauss_large.r_cut[('muc_h', 'muc_p')] = 3.5*adj_diameter
            gauss_large.r_cut[('muc_p', 'muc_p')] = 3.5*2.0*configurator.r_histone

            # Disable the excluded volume interactions from the original list
            gauss.r_cut[('muc_e', 'muc_p')] = 0.0
            gauss.r_cut[('muc_c', 'muc_p')] = 0.0
            gauss.r_cut[('muc_h', 'muc_p')] = 0.0
            gauss.r_cut[('muc_p', 'muc_p')] = 0.0

    else:

        wca = None
        ewca = None
        lje = None
        ljp = None
        eljp = None
        bmh = None

        # WCA potential for the excluded volume interactions of like-sized beads
        wca = md.pair.LJ(nlist = nlists[0])
        wca.mode='shift'

        # WCA interaction for excluded volume interaction of differently sized beads
        if configurator.nlist_n == 1:
            ewca = md.pair.ExpandedLJ(nlist = nlists[0])
            ewca.mode = 'shift'
        elif configurator.nlist_n > 1:
            ewca = md.pair.ExpandedLJ(nlist = nlists[1])
            ewca.mode = 'shift'
        else:
            print(f"ERROR: Something has gone wrong in the number of neighbor lists we are presented with, exiting!")
            sys.exit(1)

        # LJ interaction of a general form (muc_e -- muc_e)
        lje = md.pair.LJ(nlist = nlists[0])

        # LJ interaction for the attraction between muc_e and muc_p when they have the same radius
        ljp = md.pair.LJ(nlist = nlists[0])

        # Shifted LJ interaction for the attraction between muc_e and muc_p when they have different radii
        if configurator.nlist_n == 1:
            eljp = md.pair.ExpandedLJ(nlist = nlists[0])
        elif configurator.nlist_n > 1:
            eljp = md.pair.ExpandedLJ(nlist = nlists[1])

        # Born mayer huggins potential for hydrophobic interactions
        bmh = md.pair.BornMayerHuggins(nlist = nlists[0])

        # Shifted zero of a lennard jones potential so that the edge acts like a WCA or lennard jones
        # This should go to something like r == 1 at r-delta == 1, giving us the deltashift
        deltashift = 0.0
        if configurator.size_asymmetry:
            deltashift = r_histone + r_sphere - 1.0

        ###############################
        # Figure out what combination of interactions the system has
        ###############################
        # Set r_cut = 0 to disable interactions

        # muc_e <--> muc_e
        if configurator.lennard_jones_ee > 0.0:
            print(f"  muc_e <--> muc_e: LJ: epsilon={configurator.lennard_jones_ee}")
            wca.params[('muc_e', 'muc_e')]  = {'epsilon': 1.0, 'sigma': 1.0}
            lje.params[('muc_e', 'muc_e')]  = {'epsilon': configurator.lennard_jones_ee, 'sigma': 1.0}
            wca.r_cut[('muc_e', 'muc_e')]   = 0.0
            lje.r_cut[('muc_e', 'muc_e')]   = 2.5
        else:
            print(f"  muc_e <--> muc_e: WCA")
            wca.params[('muc_e', 'muc_e')]  = {'epsilon': 1.0, 'sigma': 1.0}
            lje.params[('muc_e', 'muc_e')]  = {'epsilon': 1.0, 'sigma': 1.0}
            wca.r_cut[('muc_e', 'muc_e')]   = np.power(2.0, 1/6)*1.0
            lje.r_cut[('muc_e', 'muc_e')]   = 0.0
        ewca.params[('muc_e', 'muc_e')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
        ljp.params[('muc_e', 'muc_e')]  = {'epsilon': 1.0, 'sigma': 1.0}
        eljp.params[('muc_e', 'muc_e')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
        bmh.params[('muc_e', 'muc_e')]  = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
        ewca.r_cut[('muc_e', 'muc_e')]  = 0.0
        ljp.r_cut[('muc_e', 'muc_e')]   = 0.0
        eljp.r_cut[('muc_e', 'muc_e')]  = 0.0
        bmh.r_cut[('muc_e', 'muc_e')]   = 0.0

        # muc_e <--> muc_c
        print(f"  muc_e <--> muc_c: WCA (default)")
        wca.params[('muc_e', 'muc_c')]      = {'epsilon': 1.0, 'sigma': 1.0}
        ewca.params[('muc_e', 'muc_c')]     = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
        lje.params[('muc_e', 'muc_c')]      = {'epsilon': 1.0, 'sigma': 1.0}
        ljp.params[('muc_e', 'muc_c')]      = {'epsilon': 1.0, 'sigma': 1.0}
        eljp.params[('muc_e', 'muc_c')]     = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
        bmh.params[('muc_e', 'muc_c')]      = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
        wca.r_cut[('muc_e', 'muc_c')]   = np.power(2.0, 1/6)*1.0
        ewca.r_cut[('muc_e', 'muc_c')]  = 0.0
        lje.r_cut[('muc_e', 'muc_c')]   = 0.0
        ljp.r_cut[('muc_e', 'muc_c')]   = 0.0
        eljp.r_cut[('muc_e', 'muc_c')]  = 0.0
        bmh.r_cut[('muc_e', 'muc_c')]   = 0.0

        # muc_e <--> muc_h
        print(f"  muc_e <--> muc_h: WCA (default)")
        wca.params[('muc_e', 'muc_h')]      = {'epsilon': 1.0, 'sigma': 1.0}
        ewca.params[('muc_e', 'muc_h')]     = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
        lje.params[('muc_e', 'muc_h')]      = {'epsilon': 1.0, 'sigma': 1.0}
        ljp.params[('muc_e', 'muc_h')]      = {'epsilon': 1.0, 'sigma': 1.0}
        eljp.params[('muc_e', 'muc_h')]     = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
        bmh.params[('muc_e', 'muc_h')]      = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
        wca.r_cut[('muc_e', 'muc_h')]   = np.power(2.0, 1/6)*1.0
        ewca.r_cut[('muc_e', 'muc_h')]  = 0.0
        lje.r_cut[('muc_e', 'muc_h')]   = 0.0
        ljp.r_cut[('muc_e', 'muc_h')]   = 0.0
        eljp.r_cut[('muc_e', 'muc_h')]  = 0.0
        bmh.r_cut[('muc_e', 'muc_h')]   = 0.0

        # muc_e <--> muc_p
        # This is the first time we see an adjusted radius based on the size of the muc_e and muc_p beads
        if configurator.lennard_jones > 0.0:
            if not configurator.size_asymmetry:
                print(f"  muc_e <--> muc_p: LJ: epsilon = {configurator.lennard_jones}, sigma = 1.0")
                wca.params[('muc_e', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
                ewca.params[('muc_e', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
                ljp.params[('muc_e', 'muc_p')]  = {'epsilon': configurator.lennard_jones, 'sigma': 1.0}
                eljp.params[('muc_e', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
                wca.r_cut[('muc_e', 'muc_p')]   = 0.0
                ewca.r_cut[('muc_e', 'muc_p')]  = 0.0
                ljp.r_cut[('muc_e', 'muc_p')]   = 2.5
                eljp.r_cut[('muc_e', 'muc_p')]  = 0.0
            else:
                print(f"  muc_e <--> muc_p: ELJ: epsilon = {configurator.lennard_jones}, sigma = 1.0, delta = {deltashift}")
                wca.params[('muc_e', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
                ewca.params[('muc_e', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
                ljp.params[('muc_e', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
                eljp.params[('muc_e', 'muc_p')] = {'epsilon': configurator.lennard_jones, 'sigma': 1.0, 'delta': deltashift}
                wca.r_cut[('muc_e', 'muc_p')]   = 0.0
                ewca.r_cut[('muc_e', 'muc_p')]  = 0.0
                ljp.r_cut[('muc_e', 'muc_p')]   = 0.0
                eljp.r_cut[('muc_e', 'muc_p')]  = deltashift + 2.5
        else:
            # All WCA interactions
            if not configurator.size_asymmetry:
                print(f"  muc_e <--> muc_p: WCA")
                wca.params[('muc_e', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
                ewca.params[('muc_e', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
                ljp.params[('muc_e', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
                eljp.params[('muc_e', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
                wca.r_cut[('muc_e', 'muc_p')]   = np.power(2.0, 1/6)*1.0
                ewca.r_cut[('muc_e', 'muc_p')]  = 0.0
                ljp.r_cut[('muc_e', 'muc_p')]   = 0.0
                eljp.r_cut[('muc_e', 'muc_p')]  = 0.0
            else:
                print(f"  muc_e <--> muc_p: EWCA: delta = {deltashift}")
                wca.params[('muc_e', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
                ewca.params[('muc_e', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': deltashift}
                ljp.params[('muc_e', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
                eljp.params[('muc_e', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
                wca.r_cut[('muc_e', 'muc_p')]   = 0.0
                ewca.r_cut[('muc_e', 'muc_p')]  = deltashift + np.power(2.0, 1/6)
                ljp.r_cut[('muc_e', 'muc_p')]   = 0.0
                eljp.r_cut[('muc_e', 'muc_p')]  = 0.0
        lje.params[('muc_e', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0}
        bmh.params[('muc_e', 'muc_p')] = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
        lje.r_cut[('muc_e', 'muc_p')] = 0.0
        bmh.r_cut[('muc_e', 'muc_p')] = 0.0

        # muc_c <--> muc_c
        print(f"  muc_c <--> muc_c: WCA (default)")
        wca.params[('muc_c', 'muc_c')]  = {'epsilon': 1.0, 'sigma': 1.0}
        ewca.params[('muc_c', 'muc_c')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
        lje.params[('muc_c', 'muc_c')]  = {'epsilon': 1.0, 'sigma': 1.0}
        ljp.params[('muc_c', 'muc_c')]  = {'epsilon': 1.0, 'sigma': 1.0}
        eljp.params[('muc_c', 'muc_c')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
        bmh.params[('muc_c', 'muc_c')]  = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
        wca.r_cut[('muc_c', 'muc_c')]   = np.power(2.0, 1/6)*1.0
        ewca.r_cut[('muc_c', 'muc_c')]  = 0.0
        lje.r_cut[('muc_c', 'muc_c')]   = 0.0
        ljp.r_cut[('muc_c', 'muc_c')]   = 0.0
        eljp.r_cut[('muc_c', 'muc_c')]  = 0.0
        bmh.r_cut[('muc_c', 'muc_c')]   = 0.0

        # muc_c <--> muc_h
        print(f"  muc_c <--> muc_h: WCA (default)")
        wca.params[('muc_c', 'muc_h')]  = {'epsilon': 1.0, 'sigma': 1.0}
        ewca.params[('muc_c', 'muc_h')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
        lje.params[('muc_c', 'muc_h')]  = {'epsilon': 1.0, 'sigma': 1.0}
        ljp.params[('muc_c', 'muc_h')]  = {'epsilon': 1.0, 'sigma': 1.0}
        eljp.params[('muc_c', 'muc_h')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
        bmh.params[('muc_c', 'muc_h')]  = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
        wca.r_cut[('muc_c', 'muc_h')]   = np.power(2.0, 1/6)*1.0
        ewca.r_cut[('muc_c', 'muc_h')]  = 0.0
        lje.r_cut[('muc_c', 'muc_h')]   = 0.0
        ljp.r_cut[('muc_c', 'muc_h')]   = 0.0
        eljp.r_cut[('muc_c', 'muc_h')]  = 0.0
        bmh.r_cut[('muc_c', 'muc_h')]   = 0.0

        # muc_c <--> muc_p
        if not configurator.size_asymmetry:
            print(f"  muc_c <--> muc_p: WCA")
            wca.params[('muc_c', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            ewca.params[('muc_c', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
            lje.params[('muc_c', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            ljp.params[('muc_c', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            eljp.params[('muc_c', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
            bmh.params[('muc_c', 'muc_p')]  = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            wca.r_cut[('muc_c', 'muc_p')]   = np.power(2.0, 1/6)*1.0
            ewca.r_cut[('muc_c', 'muc_p')]  = 0.0
            lje.r_cut[('muc_c', 'muc_p')]   = 0.0
            ljp.r_cut[('muc_c', 'muc_p')]   = 0.0
            eljp.r_cut[('muc_c', 'muc_p')]  = 0.0
            bmh.r_cut[('muc_c', 'muc_p')]   = 0.0
        else:
            print(f"  muc_c <--> muc_p: EWCA: delta = {deltashift}")
            wca.params[('muc_c', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            ewca.params[('muc_c', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': deltashift}
            lje.params[('muc_c', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            ljp.params[('muc_c', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            eljp.params[('muc_c', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
            bmh.params[('muc_c', 'muc_p')]  = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            wca.r_cut[('muc_c', 'muc_p')]   = 0.0
            ewca.r_cut[('muc_c', 'muc_p')]  = deltashift + np.power(2.0, 1/6)
            lje.r_cut[('muc_c', 'muc_p')]   = 0.0
            ljp.r_cut[('muc_c', 'muc_p')]   = 0.0
            eljp.r_cut[('muc_c', 'muc_p')]  = 0.0
            bmh.r_cut[('muc_c', 'muc_p')]   = 0.0

        # muc_h <--> muc_h
        if configurator.bmh < 0.0:
            print(f"  muc_h <--> muc_h: BMH: A = {configurator.bmh}, sigma = 1.0, rho = 1.0/9.0, C = 1.0, D = 1.0")
            wca.params[('muc_h', 'muc_h')] = {'epsilon': 1.0, 'sigma': 1.0}
            bmh.params[('muc_h', 'muc_h')] = {'A': configurator.bmh, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            wca.r_cut[('muc_h', 'muc_h')] = 0.0
            bmh.r_cut[('muc_h', 'muc_h')] = 2.5
        else:
            print(f"  muc_h <--> muc_h: WCA")
            wca.params[('muc_h', 'muc_h')] = {'epsilon': 1.0, 'sigma': 1.0}
            bmh.params[('muc_h', 'muc_h')] = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            wca.r_cut[('muc_h', 'muc_h')] = np.power(2.0, 1/6)*1.0
            bmh.r_cut[('muc_h', 'muc_h')] = 0.0
        ewca.params[('muc_h', 'muc_h')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
        lje.params[('muc_h', 'muc_h')]  = {'epsilon': 1.0, 'sigma': 1.0}
        ljp.params[('muc_h', 'muc_h')]  = {'epsilon': 1.0, 'sigma': 1.0}
        eljp.params[('muc_h', 'muc_h')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
        ewca.r_cut[('muc_h', 'muc_h')]  = 0.0
        lje.r_cut[('muc_h', 'muc_h')]   = 0.0
        ljp.r_cut[('muc_h', 'muc_h')]   = 0.0
        eljp.r_cut[('muc_h', 'muc_h')]  = 0.0

        # muc_h <--> muc_p
        if not configurator.size_asymmetry:
            print(f"  muc_h <--> muc_p: WCA")
            wca.params[('muc_h', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            ewca.params[('muc_h', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
            lje.params[('muc_h', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            ljp.params[('muc_h', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            eljp.params[('muc_h', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
            bmh.params[('muc_h', 'muc_p')]  = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            wca.r_cut[('muc_h', 'muc_p')]   = np.power(2.0, 1/6)*1.0
            ewca.r_cut[('muc_h', 'muc_p')]  = 0.0
            lje.r_cut[('muc_h', 'muc_p')]   = 0.0
            ljp.r_cut[('muc_h', 'muc_p')]   = 0.0
            eljp.r_cut[('muc_h', 'muc_p')]  = 0.0
            bmh.r_cut[('muc_h', 'muc_p')]   = 0.0
        else:
            print(f"  muc_h <--> muc_p: EWCA: delta = {deltashift}")
            wca.params[('muc_h', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            ewca.params[('muc_h', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': deltashift}
            lje.params[('muc_h', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            ljp.params[('muc_h', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            eljp.params[('muc_h', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
            bmh.params[('muc_h', 'muc_p')]  = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            wca.r_cut[('muc_h', 'muc_p')]   = 0.0
            ewca.r_cut[('muc_h', 'muc_p')]  = deltashift + np.power(2.0, 1/6)
            lje.r_cut[('muc_h', 'muc_p')]   = 0.0
            ljp.r_cut[('muc_h', 'muc_p')]   = 0.0
            eljp.r_cut[('muc_h', 'muc_p')]  = 0.0
            bmh.r_cut[('muc_h', 'muc_p')]   = 0.0

        # muc_c <--> muc_p
        if not configurator.size_asymmetry:
            print(f"  muc_p <--> muc_p: WCA")
            wca.params[('muc_p', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            ewca.params[('muc_p', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
            lje.params[('muc_p', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            ljp.params[('muc_p', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            eljp.params[('muc_p', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
            bmh.params[('muc_p', 'muc_p')]  = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            wca.r_cut[('muc_p', 'muc_p')]   = np.power(2.0, 1/6)*1.0
            ewca.r_cut[('muc_p', 'muc_p')]  = 0.0
            lje.r_cut[('muc_p', 'muc_p')]   = 0.0
            ljp.r_cut[('muc_p', 'muc_p')]   = 0.0
            eljp.r_cut[('muc_p', 'muc_p')]  = 0.0
            bmh.r_cut[('muc_p', 'muc_p')]   = 0.0
        else:
            print(f"  muc_p <--> muc_p: EWCA: delta = {2.0*r_histone-1.0}")
            wca.params[('muc_p', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            ewca.params[('muc_p', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 2.0*r_histone-1.0}
            lje.params[('muc_p', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            ljp.params[('muc_p', 'muc_p')]  = {'epsilon': 1.0, 'sigma': 1.0}
            eljp.params[('muc_p', 'muc_p')] = {'epsilon': 1.0, 'sigma': 1.0, 'delta': 0.0}
            bmh.params[('muc_p', 'muc_p')]  = {'A': 0.0, 'sigma': 1.0, 'rho': 1.0/9.0, 'C': 1.0, 'D': 1.0}
            wca.r_cut[('muc_p', 'muc_p')]   = 0.0
            ewca.r_cut[('muc_p', 'muc_p')]  = 2.0*r_histone - 1.0 + np.power(2.0, 1/6)
            lje.r_cut[('muc_p', 'muc_p')]   = 0.0
            ljp.r_cut[('muc_p', 'muc_p')]   = 0.0
            eljp.r_cut[('muc_p', 'muc_p')]  = 0.0
            bmh.r_cut[('muc_p', 'muc_p')]   = 0.0

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
    langevin.gamma['muc_e'] = 1.0/configurator.t_damp
    langevin.gamma['muc_c'] = 1.0/configurator.t_damp
    langevin.gamma['muc_h'] = 1.0/configurator.t_damp
    if not configurator.size_asymmetry:
        langevin.gamma['muc_p'] = 1.0/configurator.t_damp
    else:
        langevin.gamma['muc_p'] = configurator.m_histone/configurator.t_damp

    integrator.methods.append(langevin)

    # Append the forces that we need
    # Pairwise forces
    if configurator.init_type == 'create_equilibration':
        if configurator.equilibration_potential == 'gauss':
            integrator.forces.append(gauss)
        else:
            integrator.forces.append(glf)
        if configurator.nlist_n > 1:
            integrator.forces.append(gauss_large)
    elif configurator.init_type == 'production':
        integrator.forces.append(wca)
        if configurator.size_asymmetry:
            integrator.forces.append(ewca)

        if configurator.lennard_jones_ee > 0.0:
            integrator.forces.append(lje)

        if configurator.lennard_jones > 0.0:
            integrator.forces.append(ljp)
            if configurator.size_asymmetry:
                integrator.forces.append(eljp)

        if configurator.bmh < 0.0:
            integrator.forces.append(bmh)

    # FENE bond potential
    integrator.forces.append(fenewca)

    # Bending potential
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
        if configurator.equilibration_potential == 'gauss':
            logger.add(gauss, quantities=['energies', 'forces'])
        else:
            logger.add(glf, quantities=['energies', 'forces'])
        if configurator.nlist_n > 1:
            logger.add(gauss_large, quantities=['energies', 'forces'])
    logger.add(sim, quantities=['timestep', 'walltime', 'tps'])
    logger.add(thermodynamic_properties)
    
    # Display some quantities to a table while running
    output_logger = hoomd.logging.Logger(categories=['scalar', 'string'])
    status = Status(sim)
    output_logger.add(sim, quantities=['timestep', 'tps'])
    output_logger[('Status', 'etr')] = (status, 'etr', 'string')
    output_logger.add(thermodynamic_properties, quantities=['kinetic_temperature', 'pressure'])
    table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(period=configurator.nwrite_log),
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
            if configurator.equilibration_potential == 'gauss':
                gauss.params[('muc_e', 'muc_e')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
                gauss.params[('muc_e', 'muc_c')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
                gauss.params[('muc_e', 'muc_h')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
                gauss.params[('muc_e', 'muc_p')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0*adj_diameter}
                gauss.params[('muc_c', 'muc_c')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
                gauss.params[('muc_c', 'muc_h')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
                gauss.params[('muc_c', 'muc_p')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0*adj_diameter}
                gauss.params[('muc_h', 'muc_h')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0}
                gauss.params[('muc_h', 'muc_p')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0*adj_diameter}
                gauss.params[('muc_p', 'muc_p')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 2.0*configurator.r_histone}
            else:
                glf.params[('muc_e', 'muc_e')] = {'A': (iblock/nblocks*100.0), 'B': 0.0, 'r0': 1.0, 'rc': 2.0}
                glf.params[('muc_e', 'muc_c')] = {'A': (iblock/nblocks*100.0), 'B': 0.0, 'r0': 1.0, 'rc': 2.0}
                glf.params[('muc_e', 'muc_h')] = {'A': (iblock/nblocks*100.0), 'B': 0.0, 'r0': 1.0, 'rc': 2.0}
                glf.params[('muc_e', 'muc_p')] = {'A': (iblock/nblocks*100.0), 'B': 0.0, 'r0': adj_diameter, 'rc': 2.0*adj_diameter}
                glf.params[('muc_c', 'muc_c')] = {'A': (iblock/nblocks*100.0), 'B': 0.0, 'r0': 1.0, 'rc': 2.0}
                glf.params[('muc_c', 'muc_h')] = {'A': (iblock/nblocks*100.0), 'B': 0.0, 'r0': 1.0, 'rc': 2.0}
                glf.params[('muc_c', 'muc_p')] = {'A': (iblock/nblocks*100.0), 'B': 0.0, 'r0': adj_diameter, 'rc': 2.0*adj_diameter}
                glf.params[('muc_h', 'muc_h')] = {'A': (iblock/nblocks*100.0), 'B': 0.0, 'r0': 1.0, 'rc': 2.0}
                glf.params[('muc_h', 'muc_p')] = {'A': (iblock/nblocks*100.0), 'B': 0.0, 'r0': adj_diameter, 'rc': 2.0*adj_diameter}
                glf.params[('muc_p', 'muc_p')] = {'A': (iblock/nblocks*100.0), 'B': 0.0, 'r0': 2.0*r_histone, 'rc': 4.0*r_histone}

            # If we have a size asymmetry, deal with it now
            if configurator.nlist_n > 1:
                gauss_large.params[('muc_e', 'muc_p')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0*adj_diameter}
                gauss_large.params[('muc_c', 'muc_p')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0*adj_diameter}
                gauss_large.params[('muc_h', 'muc_p')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 1.0*adj_diameter}
                gauss_large.params[('muc_p', 'muc_p')] = {'epsilon': (iblock/nblocks*100.0), 'sigma': 2.0*configurator.r_histone}

            sim.run(nblock_size)

    elif configurator.init_type == 'production':
        sim.run(configurator.nsteps)

