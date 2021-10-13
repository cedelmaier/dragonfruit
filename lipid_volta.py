# Test of the grime lipid N site interaction with hoomd 300 beta9
import hoomd
import hoomd.md as md
import gsd.hoomd

import argparse
import datetime
import itertools
import os
import sys

import numpy as np

# Create a rudimentary parser to get the number of steps, the write frequency,
# and the box length
def parse_args():
    parser = argparse.ArgumentParser(prog='lipid_volta_v1.py')

    # Nsteps
    parser.add_argument('--nsteps', type=int, default=1e6,
            help='Number of steps')
    # Nwrite
    parser.add_argument('--nwrite', type=int, default=1000,
            help='Number of steps between output')

    # Various simulation quantities
    parser.add_argument('--lbox', type=float, default=100.0,
            help='Box extent (sigma, sim units)')
    parser.add_argument('--ngrid', type=int, default=100,
            help='Box grid points')

    # Lipid parameters
    parser.add_argument('--mdopc', type=float, default=734.039,
            help='Lipid mass (m, sim units)')
    parser.add_argument('--lgamma', type=float, default=0.872,
            help='Lipid drag coefficient (m/tau, sim units)')
    parser.add_argument('--lA', type=float, default=187.5,
            help='Lipid repulsion coeff. A (kT/sigma, sim units)')
    parser.add_argument('--lB', type=float, default=7.0/3.0,
            help='Lipid attraction coeff. B (kT/sigma, sim units)')
    parser.add_argument('--lkbond', type=float, default=2812.5,
            help='Lipid bond k (kT/sigma^2, sim units)')
    parser.add_argument('--lkbend', type=float, default=2.0,
            help='Lipid bend k (kT/rad^2, sim units)')

    # Add in the Ah parameters
    parser.add_argument('--nah', type=int, default=0,
            help='Number of AH domains')
    parser.add_argument('--mah', type=float, default=2535.3864,
            help='AH mass (m, sim units)')
    parser.add_argument('--ahgamma', type=float, default=0.872,
            help='Lipid drag coefficient (m/tau, sim units)')
    parser.add_argument('--ahA', type=float, default=187.5,
            help='AH domain repulsion coeff. A (kT/sigma, sim units)')
    parser.add_argument('--ahB', type=float, default=7.0/3.0,
            help='AH domain attraction coeff. B (kT/sigma, sim units)')
    parser.add_argument('--ahkbond', type=float, default=2812.5,
            help='AH bond k (kT/sigma^2, sim units)')
    parser.add_argument('--ahkbend', type=float, default=2.0,
            help='AH bend k (kT/rad^2, sim units)')

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

    # Simulation parameters
    kT = 1.0
    bead_size = 1.0
    deltatau = 0.05
    nsteps = np.int32(opts.nsteps)
    nwrite = np.int32(opts.nwrite)
    ngrid = np.int32(opts.ngrid)
    lbox = np.float32(opts.lbox)

    # Lipid parameters
    nbeads = 4
    mdopc = np.float32(opts.mdopc)
    lgamma = np.float32(opts.lgamma)
    lA = np.float32(opts.lA)
    lB = np.float32(opts.lB)
    lkbond = np.float32(opts.lkbond)
    lkbend = np.float32(opts.lkbend)

    # AH-domain parameters
    ahgamma = np.float32(opts.lgamma)
    ahmass = np.float32(opts.mah)
    ahA = np.float32(opts.ahA)
    ahB = np.float32(opts.ahB)
    ahkbond = np.float32(opts.ahkbond)
    ahkbend = np.float32(opts.ahkbend)

    # Create the hoomd context on the GPU
    gpu = hoomd.device.GPU(notice_level = 2)
    sim = hoomd.Simulation(gpu, seed=1)

    # Print the system information
    print(f"\n--------")
    print(f"System information")
    print(f"Simulation time (tau)   = {deltatau*nsteps}")
    print(f"kBT                     = {kT}")
    
    ###############################################################################
    # Setup the system
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
    print(f"\n--------")
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
    ah_nwidth = 2 # 2 beads thick
    ah_nlength = 5 # 10 beads long
    ah_diameter = 20./7.5 # Diameter of ah in sigma
    ah_bead_size = ah_diameter/ah_nwidth
    ahmass_per_bead = ahmass / (ah_nwidth * ah_nlength)
    ahr0 = ah_bead_size
    ahrc = 2.0*ah_bead_size
    ahrbond = ahr0

    # Print out AH information
    print(f"\n--------")
    print(f"AH information")
    print(f"AH bead size (sigma, sim units)     = {ah_bead_size}")
    print(f"AH mass (amu)                       = {ahmass}")
    print(f"  AH mass per bead                  = {ahmass_per_bead}")
    print(f"Gamma (m/tau)                       = {ahgamma}")
    print(f"AH A (kBT / sigma)                  = {ahA}")
    print(f"AH B (kBT / sigma)                  = {ahB}")
    print(f"AH R (sigma)                        = {ahr0}")
    print(f"AH RC (sigma)                       = {ahrc}")
    print(f"AH r (sigma)                        = {ahrbond}")
    print(f"AH k (kBT/sigma^2)                  = {ahkbond}")
    print(f"AH k (kBT/rad^2)                    = {ahkbend}")

    ###############################################################################
    # Creation of the system
    ###############################################################################
    # What is the number of replications we need in the box?
    box_extent = lbox
    linear_extent = box_extent/ngrid
    fudge_size = 0.2 # Prevent particles on the boundaries because reasons
    
    # Create the 2 bilayers, so double the size of everything
    snap = hoomd.Snapshot(gpu.communicator)
    snap.configuration.box = [linear_extent, linear_extent, 2.0*((bead_size/2.0) + (nbeads-1)*bead_size)+fudge_size, 0, 0, 0]
    
    # Change how we configure the system. Numbers run from 0 to nbeads, then reverse the z-direction and do the mirror image
    # plane. This should be more robust than previous versions to handle variable number of beads
    
    # Create the number of bead types we need (including AH domains)
    snap.particles.N = 2*nbeads
    snap.particles.types = ['H', 'I', 'T', 'AHH', 'AHT']
    
    # The number of bonds we have is determined by nbeads, as is the number of angle potentials
    snap.bonds.N = 2*(nbeads - 1)
    snap.angles.N = 2*(nbeads - 2)
    snap.bonds.types = ['lipidbond', 'ahbond']
    snap.bonds.typeid[:] = 0
    snap.angles.types = ['lipidbend', 'ahbend', 'ahbend90']
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
    
    # Create the ah filaments themselves after a box update
    ah_start_nx = snap.particles.N
    snap.particles.N = snap.particles.N + (ah_nlength * ah_nwidth)
    # Assign the location of the AH filaments
    # First row of filaments
    ah_start_x = np.array([0.0, 0.0, box_extent/4.0])
    for idx in range(ah_start_nx, snap.particles.N - ah_nlength):
        snap.particles.position[idx] = ah_start_x + np.array([0, (ah_bead_size/2.0)*(idx - ah_start_nx), 0])
        snap.particles.typeid[idx] = 3 # AHH
        snap.particles.mass[idx] = ahmass_per_bead
    # Second row of filaments
    for idx in range(ah_start_nx + ah_nlength, snap.particles.N):
        snap.particles.position[idx] = ah_start_x + np.array([0, (ah_bead_size/2.0)*(idx - ah_start_nx - ah_nlength), -(ah_bead_size/2.0)])
        snap.particles.typeid[idx] = 4 # AHT
        snap.particles.mass[idx] = ahmass_per_bead
    
    # Print the information
    for idx in range(ah_start_nx, snap.particles.N):
        print("Particle ID {}".format(idx))
        print("  position: {}".format(snap.particles.position[idx]))
    
    # Now set up the bond information
    n_newbonds = 1 + 3*(ah_nlength-1)
    ah_start_bdx = snap.bonds.N
    bdx = ah_start_bdx
    snap.bonds.N = snap.bonds.N + n_newbonds
    for idx in range(ah_start_nx, snap.particles.N - ah_nlength - 1):
        snap.bonds.typeid[bdx] = 1
        snap.bonds.group[bdx] = [idx, idx+1]
        bdx += 1
    for idx in range(ah_start_nx, snap.particles.N - ah_nlength):
        snap.bonds.typeid[bdx] = 1
        snap.bonds.group[bdx] = [idx, idx+ah_nlength]
        bdx += 1
    for idx in range(ah_start_nx, snap.particles.N - ah_nlength - 1):
        snap.bonds.typeid[bdx] = 1
        snap.bonds.group[bdx] = [idx+ah_nlength, idx+ah_nlength+1]
        bdx += 1
    
    for bdx in range(ah_start_bdx, snap.bonds.N):
        print("Bond ID: {}".format(bdx))
        print("  group: {}".format(snap.bonds.group[bdx]))
    
    # Set up the angle information
    n_newangles = 2*(ah_nlength - 2)
    ah_start_adx = snap.angles.N
    adx = ah_start_adx
    snap.angles.N = snap.angles.N + n_newangles
    for idx in range(ah_start_nx + 1, ah_start_nx + 4):
        print("Angle: {}".format(adx))
        snap.angles.typeid[adx] = 1
        print("  type: {}".format(snap.angles.typeid[adx]))
        snap.angles.group[adx] = [idx-1,idx,idx+1]
        print("  group1: {}".format(snap.angles.group[adx]))
        adx += 1
        snap.angles.group[adx]= [idx-1+ah_nlength,idx+ah_nlength,idx+1+ah_nlength]
        print("  group2: {}".format(snap.angles.group[adx]))
        adx += 1
    
    # Add another new angle for every group of 4 that we have
    n_newangles = 4 * (ah_nlength - 1)
    ah_start_adx = snap.angles.N
    adx = ah_start_adx
    snap.angles.N = snap.angles.N + n_newangles
    for idx in range(ah_start_nx, ah_start_nx + ah_nlength-1):
        snap.angles.typeid[adx+0] = 2
        snap.angles.typeid[adx+1] = 2
        snap.angles.typeid[adx+2] = 2
        snap.angles.typeid[adx+3] = 2
    
        # 0 1 6
        # 1 6 5
        # 6 5 0
        # 5 0 1
        snap.angles.group[adx+0] = [idx, idx+1, idx+1+ah_nlength]
        snap.angles.group[adx+1] = [idx+1, idx+1+ah_nlength, idx+ah_nlength]
        snap.angles.group[adx+2] = [idx+1+ah_nlength, idx+ah_nlength, idx]
        snap.angles.group[adx+3] = [idx+ah_nlength, idx, idx+1]
    
        for i in range(0,4):
            print("Angle: {}".format(adx+i))
            print("  type: {}".format(snap.angles.typeid[adx+i]))
            print("  group: {}".format(snap.angles.group[adx+i]))
    
        adx = adx + 4
    
    ###############################################################################
    # Create the system
    ###############################################################################
    sim.create_state_from_snapshot(snap)
    
    # Create the pair potential
    glf = md.pair.GrimeLipid(nlist=md.nlist.Cell())
    
    # Setting the particles is interesting. First off, assign zero to every combination
    # in the system
    lipid_types = ['H', 'I', 'T']
    ah_types = ['AHH', 'AHT']
    
    ###############################
    # Membrane-membrane interactions
    ###############################
    # Assign all of the coefficients between particles
    #Head to other interactions
    glf.params[('H', 'H')] = {'A': lAhh, 'B': 0.0, 'r0': lrhh, 'rc': lrhh} # Head-Head
    glf.params[('H', 'I')] = {'A': lA, 'B': 0.0, 'r0': lr0, 'rc': lrc} # Head-Interface
    glf.params[('H', 'T')] = {'A': lA, 'B': 0.0, 'r0': lr0, 'rc': lrc} # Head-Tail
    
    #Interface to others interactions
    glf.params[('I', 'I')] = {'A': lA, 'B': lB, 'r0': lr0, 'rc': lrc} # Interface to interface is sticky
    glf.params[('I', 'T')] = {'A': lA, 'B': 0.0, 'r0': lr0, 'rc': lrc} # Interface to tail is not sticky
    
    #Tail to other interactions
    glf.params[('T', 'T')] = {'A': lA, 'B': lB, 'r0': lr0, 'rc': lrc} # Tail to tail is sticky
    
    # Set the cutoff distance accordingly
    # FIXME: Check this, as now we can have different exclusions
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
    # system.  Self interactions for the ahs is zero, this should be taken care of in the
    # bond and bend potentials
    glf.params[('AHH', 'AHH')] = {'A': ahA, 'B': 0.0, 'r0': ahr0, 'rc': ahrc}
    glf.params[('AHH', 'AHT')] = {'A': ahA, 'B': 0.0, 'r0': ahr0, 'rc': ahrc}
    glf.params[('AHT', 'AHT')] = {'A': ahA, 'B': 0.0, 'r0': ahr0, 'rc': ahrc}
    glf.r_cut[('AHH', 'AHH')] = ahr0
    glf.r_cut[('AHH', 'AHT')] = ahr0
    glf.r_cut[('AHT', 'AHT')] = ahr0
    
    ###############################
    # Membrane-AH interactions
    ###############################
    # Set all the potentials to be purely repulsive between AH and the membrane
    # beads
    for x in itertools.product(lipid_types, ah_types):
        glf.params[x] = {'A': ahA, 'B': 0.0, 'r0': 0.5*(ahr0+lr0), 'rc': ahrc}
        glf.r_cut[x] = 0.5*(ahr0+lr0)
    
    # Set membrane-head ah-head interaction
    glf.params[('H', 'AHH')] = {'A': ahA, 'B': ahB, 'r0': 0.5*(ahr0+lr0), 'rc': 1.0*(ahr0+lr0)}
    glf.r_cut[('H', 'AHH')] = 1.0*(ahr0+lr0)
    # Now set only the tail-tail interactions for ahs and the membrane
    glf.params[('T', 'AHT')] = {'A': ahA, 'B': ahB, 'r0': 0.5*(ahr0+lr0), 'rc': 1.0*(ahr0+lr0)}
    glf.r_cut[('T', 'AHT')] = 1.0*(ahr0+lr0)
    
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
    angleharmonic.params['ahbend'] = dict(k = ahkbend, t0 = np.pi)
    angleharmonic.params['ahbend90'] = dict(k = ahkbend, t0 = np.pi/2.0)
    
    # Set up the integrator for the system
    #nph = hoomd.md.methods.NPH(filter = hoomd.filter.All(),
    #                           tauS = pdamp_per_tau,
    #                           S = 0.0,
    #                           box_dof = [True, True, False, False, False, False],
    #                           couple = "xy")
    langevin = hoomd.md.methods.Langevin(hoomd.filter.All(), kT = kT)
    langevin.gamma['H'] = lgamma
    langevin.gamma['I'] = lgamma
    langevin.gamma['T'] = lgamma
    langevin.gamma['AHH'] = ahgamma
    langevin.gamma['AHT'] = ahgamma
    
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
    
    ###############################################################################
    # Write out any interesting information
    ###############################################################################
    nlipids = ngrid * ngrid * 2
    nlipids_per_leafleat = nlipids / 2
    area_per_lipid = (sim.state.box.Lx * sim.state.box.Ly) / (nlipids_per_leafleat)
    print(f"System paramters:")
    print(f"  Number of particles:                  {sim.state.N_particles}")
    print(f"  Number of lipid patches:              {nlipids}")
    print(f"  Number of lipid patches per leafleat: {nlipids_per_leafleat}")
    print(f"  Area per lipid:                       {area_per_lipid}")
    print(f"  Box size:                             {sim.state.box}")
    
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
    fname_gsd = 'traj_volta_v1_lbox{:4.2f}_ngrid{:d}_lA{:4.2f}_lB{:4.2f}_ahA{:4.2f}_ahB{:4.2f}.gsd'.format(lbox, ngrid, lA, lB, ahA, ahB)
    gsd_writer = hoomd.write.GSD(filename = fname_gsd,
                                 trigger = hoomd.trigger.Periodic(nwrite),
                                 mode = 'wb',
                                 filter = hoomd.filter.All(),
                                 log = logger)
    sim.operations.writers.append(gsd_writer)
    
    sim.run(nsteps)
    
