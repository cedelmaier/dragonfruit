# Test of the grime lipid N site interaction with hoomd 300 beta9
import hoomd
import hoomd.md as md
import gsd.hoomd

import argparse
import datetime
import itertools
import os
import sys
import yaml

import numpy as np

from Configurator import Configurator

# Create a rudimentary parser to get the number of steps, the write frequency,
# and the box length
def parse_args():
    parser = argparse.ArgumentParser(prog='lipid_volta.py')

    parser.add_argument('--default_file', type=str, default='lipid.default.yaml',
            help='Default configuration file')

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

    configurator = Configurator(opts)

    # Simulation parameters
    kT = np.float32(configurator.default_yaml['simulation']['kT'])
    bead_size = np.float32(configurator.default_yaml['simulation']['bead_size'])
    deltatau = np.float32(configurator.default_yaml['simulation']['deltatau'])
    nsteps = np.int32(np.float32(configurator.default_yaml['simulation']['nsteps']))
    nwrite = np.int32(np.float32(configurator.default_yaml['simulation']['nwrite']))
    ngrid = np.int32(np.float32(configurator.default_yaml['simulation']['ngrid']))
    lbox = np.float32(configurator.default_yaml['simulation']['lbox'])
    nseed = np.int32(np.float32(configurator.default_yaml['simulation']['seed']))
    compute_mode = configurator.default_yaml['simulation']['mode']
    trajectory_file = configurator.default_yaml['simulation']['trajectory_file']

    # Lipid parameters
    nbeads = np.int32(np.float32(configurator.default_yaml['membrane']['nbeads']))
    mdopc = np.float32(configurator.default_yaml['membrane']['mass'])
    lgamma = np.float32(configurator.default_yaml['membrane']['gamma'])
    lA = np.float32(configurator.default_yaml['membrane']['A'])
    lB = np.float32(configurator.default_yaml['membrane']['B'])
    lkbond = np.float32(configurator.default_yaml['membrane']['kbond'])
    lkbend = np.float32(configurator.default_yaml['membrane']['kbend'])

    # AH-domain parameters
    n_ah = np.int32(np.float32(configurator.default_yaml['ah_domain']['n_ah']))
    polymer_type = configurator.default_yaml['ah_domain']['polymer_type']
    ahgamma = np.float32(configurator.default_yaml['ah_domain']['gamma'])
    ahnbeads = np.int32(np.float32(configurator.default_yaml['ah_domain']['nbeads']))
    ahnrepeat = np.int32(np.float32(configurator.default_yaml['ah_domain']['nrepeat']))
    ahmass = np.float32(configurator.default_yaml['ah_domain']['mass'])
    ahA = np.float32(configurator.default_yaml['ah_domain']['A'])
    ahBsurface = np.float32(configurator.default_yaml['ah_domain']['B_surface'])
    ahBintermediate = np.float32(configurator.default_yaml['ah_domain']['B_intermediate'])
    ahBdeep = np.float32(configurator.default_yaml['ah_domain']['B_deep'])
    ahBself = np.float32(configurator.default_yaml['ah_domain']['B_self'])
    ahkbond = np.float32(configurator.default_yaml['ah_domain']['kbond'])

    # Create the hoomd context on the GPU
    if compute_mode == 'cpu':
        cpu = hoomd.device.CPU(notice_level = 2)
        sim = hoomd.Simulation(cpu, seed = nseed)
    elif compute_mode == 'gpu':
        gpu = hoomd.device.GPU(notice_level = 2)
        sim = hoomd.Simulation(gpu, seed = nseed)

    # Print the system information
    print(f"--------")
    print(f"System information")
    print(f"Simulation time (tau)   = {deltatau*nsteps}")
    print(f"kBT                     = {kT}")
    print(f"seed                    = {nseed}")
    print(f"box size:               = {lbox}")
    
    ###############################################################################
    # Set up the system
    ###############################################################################
    # Head beads are slightly smaller than the normal beads
    lr0 = bead_size
    lrc = 2.0 * lr0
   
    # Special choices for the head bead size
    lrhh = 0.75 * lr0
    lAhh = 0.75 * lA
    
    # Bond distance and length, set to the distance between molecules
    lrbond = lr0
    
    # Lipid mass in g/mol 
    lipidmass = mdopc
    lipidmass_per_bead = lipidmass / nbeads
    
    # Print out the relevant information about the lipids
    print(f"--------")
    print(f"Lipid information")
    print(f"Lipid mass (amu)        = {lipidmass}")
    print(f"  Lipid mass per bead   = {lipidmass_per_bead}")
    print(f"Gamma (m/tau)           = {lgamma}")
    print(f"A (kBT)                 = {lA}")
    print(f"B (kBT)                 = {lB}")
    print(f"R (sigma)               = {lr0}")
    print(f"RC (sigma)              = {lrc}")
    print(f"(Ahh (kBT))             = {lAhh}")
    print(f"(Rhh (sigma))           = {lrhh}")
    print(f"Bond r (sigma)          = {lrbond}")
    print(f"Bond k (kBT/sigma^2)    = {lkbond}")
    print(f"Bend k (kBT/rad^2)      = {lkbend}")
    
    # Create the ah single dimer itself (for now)
    ah_nlength = 24 # 24 beads long
    ah_diameter = 0.75 # 5.625 angstroms diameter for each bead
    ah_bead_size = ah_diameter
    ahmass_per_bead = ahmass / (ah_nlength)
    ahr0 = ah_bead_size
    ahrc = 2.0*ah_bead_size
    ahrbond = ahr0

    # Print out AH information
    print(f"--------")
    print(f"AH information (copolymer model)")
    print(f"AH number                           = {n_ah}")
    print(f"AH bead size (sigma)                = {ah_bead_size}")
    print(f"AH mass (amu)                       = {ahmass}")
    print(f"  AH mass per bead                  = {ahmass_per_bead}")
    print(f"Gamma (m/tau)                       = {ahgamma}")
    print(f"AH A (kBT / sigma)                  = {ahA}")
    print(f"AH B terms (kBT / sigma)")
    print(f"  B surface                         = {ahBsurface}")
    print(f"  B intermediate                    = {ahBintermediate}")
    print(f"  B deep                            = {ahBdeep}")
    print(f"  B self                            = {ahBself}")
    print(f"AH R (sigma)                        = {ahr0}")
    print(f"AH RC (sigma)                       = {ahrc}")
    print(f"AH r (sigma)                        = {ahrbond}")
    print(f"AH k (kBT/sigma^2)                  = {ahkbond}")

    # Write out any interesting derived information
    nlipids = ngrid * ngrid * 2
    nlipids_per_leafleat = nlipids / 2
    area_per_lipid = (lbox * lbox) / (nlipids_per_leafleat)
    print(f"--------")
    print(f"Derived information")
    print(f"  Number of lipid patches:              {nlipids}")
    print(f"  Number of lipid patches per leafleat: {nlipids_per_leafleat}")
    print(f"  Area per lipid:                       {area_per_lipid}")
    

    ###############################################################################
    # Creation of the system
    ###############################################################################
    # What is the number of replications we need in the box?
    box_extent = lbox
    linear_extent = box_extent/ngrid
    fudge_size = 0.2 # Prevent particles on the boundaries because reasons
    
    # Create the 2 bilayers, so double the size of everything
    if compute_mode == 'cpu':
        snap = hoomd.Snapshot(cpu.communicator)
    elif compute_mode == 'gpu':
        snap = hoomd.Snapshot(gpu.communicator)

    snap.configuration.box = [linear_extent, linear_extent, 2.0*((bead_size/2.0) + (nbeads-1)*bead_size)+fudge_size, 0, 0, 0]
    
    # Change how we configure the system. Numbers run from 0 to nbeads, then reverse the z-direction and do the mirror image
    # plane. This should be more robust than previous versions to handle variable number of beads
    
    # Create the number of bead types we need (including AH domains)
    snap.particles.N = 2*nbeads
    snap.particles.types = ['H', 'I', 'T', 'AH1', 'AH2']
    
    # The number of bonds we have is determined by nbeads, as is the number of angle potentials
    snap.bonds.N = 2*(nbeads - 1)
    snap.angles.N = 2*(nbeads - 2)
    snap.bonds.types = ['lipidbond', 'ahbond']
    snap.bonds.typeid[:] = 0
    snap.angles.types = ['lipidbend']
    snap.angles.typeid[:] = 0
    
    # We have two leaflets, assign each with a zleaf loop
    # Start idx at 0, this will increment for every change!
    idx = 0
    for zidx in [1, -1]:
        for i in range(nbeads):
            snap.particles.position[idx] = [0.0, 0.0, zidx*( (bead_size/2.0) + (nbeads-1-i)*bead_size )]
    
            # Set the typeid in an ugly if else block, but easiest way to control the double Tail portion
            if i == 0:
                snap.particles.typeid[idx] = 0
            elif i == 1:
                snap.particles.typeid[idx] = 1
            else:
                snap.particles.typeid[idx] = 2
    
            # Set the mass for the beads to be the same
            snap.particles.mass[idx] = lipidmass_per_bead
    
            idx += 1
    
    # Sadly, hardcoding the bond types for each number of beads seems best
    if nbeads == 3:
        snap.bonds.group[:] = [[0, 1], [1, 2], [3, 4], [4, 5]]
        snap.angles.group[:] = [[0, 1, 2], [3, 4, 5]]
    elif nbeads == 4:
        snap.bonds.group[:] = [[0, 1], [1, 2], [2, 3], [4, 5], [5, 6], [6, 7]]
        snap.angles.group[:] = [[0, 1, 2], [1, 2, 3], [4, 5, 6], [5, 6, 7]]
    else:
        print("Only used for 3,4-bead models, exiting", file=sys.stderr)
    
    # Replicate the snapshot on a lattice
    snap.replicate(ngrid, ngrid, 1)
    
    # Create the box size of the system to be something reasonable
    snap.configuration.box = [box_extent, box_extent, box_extent, 0, 0, 0]

    # XXX: Move the configuration into a separate script
    if n_ah != 0:
        # Create the ah filaments themselves after a box update
        # This is for the copolymer model, some linear number of bonds
        ah_start_nx = snap.particles.N
        snap.particles.N = snap.particles.N + ah_nlength
        # Assign the locatin of the AH filament
        ah_start_x = np.array([0.0, 0.0, box_extent/4.0])
        ndx = 0
        ah_code = False
        for idx in range(ah_start_nx, snap.particles.N):
            if ndx % 3 == 0:
                ah_code = not ah_code
            snap.particles.position[idx] = ah_start_x + np.array([0, (ah_bead_size/2.0)*(idx - ah_start_nx), 0])
            if ah_code:
                snap.particles.typeid[idx] = 3 # AH1
            else:
                snap.particles.typeid[idx] = 4 # AH2
            snap.particles.mass[idx] = ahmass_per_bead

            ndx += 1

        #for idx in range(ah_start_nx, snap.particles.N):
        #    print("ID: {}".format(idx))
        #    print("  typeid: {}, position: {}".format(snap.particles.typeid[idx], snap.particles.position[idx]))

        # Now set up the bond information
        n_newbonds = ah_nlength - 1
        ah_start_bdx = snap.bonds.N
        bdx = ah_start_bdx
        snap.bonds.N = snap.bonds.N + n_newbonds
        for idx in range(ah_start_nx, snap.particles.N - 1):
            snap.bonds.typeid[bdx] = 1
            snap.bonds.group[bdx] = [idx, idx+1]
            bdx += 1
        
        #for bdx in range(ah_start_bdx, snap.bonds.N):
        #    print("Bond ID: {}".format(bdx))
        #    print("  group: {}".format(snap.bonds.group[bdx]))
    
    ###############################################################################
    # Create the system
    ###############################################################################
    sim.create_state_from_snapshot(snap)
    
    # Create the pair potential
    glf = md.pair.GrimeLipid(nlist=md.nlist.Cell())
    
    # Setting the particles is interesting. First off, assign zero to every combination
    # in the system
    lipid_types = ['H', 'I', 'T']
    ah_types = ['AH1', 'AH2']
    
    ###############################
    # Membrane-membrane interactions
    ###############################
    # Assign all of the coefficients between particles
    #Head to other interactions
    glf.params[('H', 'H')] = {'A': lAhh, 'B': 0.0, 'r0': lrhh, 'rc': 2.0*lrhh} # Head-Head
    glf.params[('H', 'I')] = {'A': lA, 'B': 0.0, 'r0': lr0, 'rc': lrc} # Head-Interface
    glf.params[('H', 'T')] = {'A': lA, 'B': 0.0, 'r0': lr0, 'rc': lrc} # Head-Tail
    
    #Interface to others interactions
    glf.params[('I', 'I')] = {'A': lA, 'B': lB, 'r0': lr0, 'rc': lrc} # Interface to interface is sticky
    glf.params[('I', 'T')] = {'A': lA, 'B': 0.0, 'r0': lr0, 'rc': lrc} # Interface to tail is not sticky
    
    #Tail to other interactions
    glf.params[('T', 'T')] = {'A': lA, 'B': lB, 'r0': lr0, 'rc': lrc} # Tail to tail is sticky
    
    # Set the cutoff distance accordingly
    glf.r_cut[('H', 'H')] = lrhh
    glf.r_cut[('H', 'I')] = lr0
    glf.r_cut[('H', 'T')] = lr0
    glf.r_cut[('I', 'I')] = lrc
    glf.r_cut[('I', 'T')] = lr0
    glf.r_cut[('T', 'T')] = lrc
    
    ###############################
    # AH-AH interactions
    ###############################
    # We also have to assign interactions between the ah pieces and every other body in the
    # system.  Self interactions for the ahs is not always zero!
    glf.params[('AH1', 'AH1')] = {'A': ahA, 'B': ahBself, 'r0': ahr0, 'rc': ahrc}
    glf.params[('AH1', 'AH2')] = {'A': ahA, 'B': 0.0, 'r0': ahr0, 'rc': ahrc}
    glf.params[('AH2', 'AH2')] = {'A': ahA, 'B': ahBself, 'r0': ahr0, 'rc': ahrc}
    glf.r_cut[('AH1', 'AH1')] = ahrc
    glf.r_cut[('AH1', 'AH2')] = ahr0
    glf.r_cut[('AH2', 'AH2')] = ahrc
    
    ###############################
    # Membrane-AH interactions
    ###############################
    # Initialize interactions to be purely repulsive to start with, use the average distance
    # as we have particles of different size interacting
    for x in itertools.product(lipid_types, ah_types):
        glf.params[x] = {'A': ahA, 'B': 0.0, 'r0': 0.5*(ahr0+lr0), 'rc': 1.0*(ahr0+lr0)}
        glf.r_cut[x] = 0.5*(ahr0+lr0)
    
    # Set membrane-head ah-head interaction
    glf.params[('H', 'AH1')] = {'A': ahA, 'B': ahBsurface, 'r0': 0.5*(ahr0+lr0), 'rc': 1.0*(ahr0+lr0)}
    glf.r_cut[('H', 'AH1')] = 1.0*(ahr0+lr0)
    # Set the intermediate range interaction for AH-second
    glf.params[('I', 'AH2')] = {'A': ahA, 'B': ahBintermediate, 'r0': 0.5*(ahr0+lr0), 'rc': 1.0*(ahr0+lr0)}
    glf.r_cut[('I', 'AH2')] = 1.0*(ahr0+lr0)
    # Now set only the tail-tail interactions for ahs and the membrane
    glf.params[('T', 'AH2')] = {'A': ahA, 'B': ahBdeep, 'r0': 0.5*(ahr0+lr0), 'rc': 1.0*(ahr0+lr0)}
    glf.r_cut[('T', 'AH2')] = 1.0*(ahr0+lr0)
    
    ###############################
    # Bonded and angle interactions
    ###############################
    # Assign bonded interaction strengths
    harmonic = md.bond.Harmonic()
    harmonic.params['lipidbond'] = dict(k = lkbond, r0 = lrbond)
    harmonic.params['ahbond'] = dict(k = ahkbond, r0 = ahrbond)
    
    # Assign angle interaction strengths
    angleharmonic = md.angle.Harmonic()
    angleharmonic.params['lipidbend'] = dict(k = lkbend, t0 = np.pi)
    
    # Set up the integrator for the system
    #nph = hoomd.md.methods.NPH(filter = hoomd.filter.All(),
    #                           tauS = pdamp,
    #                           S = 0.0,
    #                           box_dof = [True, True, False, False, False, False],
    #                           couple = "xy")
    langevin = hoomd.md.methods.Langevin(hoomd.filter.All(),
                                         kT = kT)
    langevin.gamma['H'] = lgamma
    langevin.gamma['I'] = lgamma
    langevin.gamma['T'] = lgamma
    langevin.gamma['AH1'] = ahgamma
    langevin.gamma['AH2'] = ahgamma
    
    integrator = md.Integrator(dt = deltatau)
    #integrator.methods.append(nph)
    integrator.methods.append(langevin)
    
    integrator.forces.append(glf)
    integrator.forces.append(harmonic)
    integrator.forces.append(angleharmonic)
    
    # Add to the simulation
    sim.operations.integrator = integrator
    
    # Run a simulation to get variables to work out, also thermalize the momenta
    sim.run(0)
    sim.state.thermalize_particle_momenta(hoomd.filter.All(), kT)
    #nph.thermalize_barostat_dof()
    
    ###############################################################################
    # Print information for the main program
    ###############################################################################
    # Keep track of the thermodynamic information
    thermodynamic_properties = md.compute.ThermodynamicQuantities(
            filter = hoomd.filter.All())
    sim.operations.computes.append(thermodynamic_properties)
    
    # Set up the logging of important quantities
    logger = hoomd.logging.Logger()
    logger.add(glf, quantities=['energies', 'forces'])
    logger.add(sim, quantities=['timestep', 'walltime', 'tps'])
    logger.add(thermodynamic_properties)
    
    # Display some quantities to a table while running
    output_logger = hoomd.logging.Logger(categories=['scalar', 'string'])
    status = Status(sim)
    output_logger.add(sim, quantities=['timestep', 'tps'])
    output_logger[('Status', 'etr')] = (status, 'etr', 'string')
    output_logger.add(thermodynamic_properties, quantities=['kinetic_temperature', 'pressure'])
    table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(period=nwrite),
                              logger=output_logger)
    sim.operations.writers.append(table)
    
    # Set up writing out to a GSD file for trajectories
    gsd_writer = hoomd.write.GSD(filename = trajectory_file,
                                 trigger = hoomd.trigger.Periodic(nwrite),
                                 mode = 'wb',
                                 filter = hoomd.filter.All(),
                                 log = logger)
    sim.operations.writers.append(gsd_writer)
    
    print(f"--------")
    print(f"Running...")
    sim.run(nsteps)
    
