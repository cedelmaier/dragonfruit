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

# Magic to get the library directory properly
sys.path.append(os.path.join(os.path.dirname(__file__), 'septins'))
from septin_seed import SeptinSeed

# Create a rudimentary parser to get the number of steps, the write frequency,
# and the box length
def parse_args():
    parser = argparse.ArgumentParser(prog='LipidVolta.py')

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
    configurator = SeptinSeed(opts.workdir, opts)

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

    # Check how we are reading in information
    if configurator.init_type == 'read_gsd':
        ftraj = gsd.hoomd.open(configurator.init_filename, 'rb')
        gsd_snap = ftraj[-1]
        ftraj.close()
        # Convert to a hoomd snapshot
        if configurator.compute_mode == 'cpu':
            snap = snap.from_gsd_snapshot(gsd_snap, cpu.communicator)
        elif configurator.compute_mode == 'gpu':
            snap = snap.from_gsd_snapshot(gsd_snap, gpu.communicator)
    elif configurator.init_type == 'all':
        fudge_size = 0.2 # Prevent particles on the boundaries because reasons
    else:
        print(f"Configurator init {configurator.init_type} not yet available, exiting!")
        sys.exit(1)

    # Initialize the membrane and AH domains with information from the snapshot (if available)
    configurator.Configure(snap)

    # Print information
    configurator.PrintInformation(snap)

    ###############################################################################
    # Create the system
    ###############################################################################
    sim.create_state_from_snapshot(snap)

    # Create the pair potential
    cell = hoomd.md.nlist.Cell(buffer=0.4)
    glf = md.pair.GrimeLipid(nlist=cell)
    
    # Setting the particles is interesting. First off, assign zero to every combination
    # in the system
    ss_types = ['SS1']
    lipid_types = ['H', 'I', 'T']
    ah_types = ['AH1', 'AH2']

    if configurator.lipids:
        ###############################
        # Membrane-membrane interactions
        ###############################
        # Assign all of the coefficients between particles
        #Head to other interactions
        glf.params[('H', 'H')] = {'A': configurator.lipids.Ahh, 'B': 0.0, 'r0': configurator.lipids.rhh, 'rc': 2.0*configurator.lipids.rhh} # Head-Head
        glf.params[('H', 'I')] = {'A': configurator.lipids.A, 'B': 0.0, 'r0': configurator.lipids.r0, 'rc': configurator.lipids.rc} # Head-Interface
        glf.params[('H', 'T')] = {'A': configurator.lipids.A, 'B': 0.0, 'r0': configurator.lipids.r0, 'rc': configurator.lipids.rc} # Head-Tail
        
        #Interface to others interactions
        glf.params[('I', 'I')] = {'A': configurator.lipids.A, 'B': configurator.lipids.B, 'r0': configurator.lipids.r0, 'rc': configurator.lipids.rc} # Interface to interface is sticky
        glf.params[('I', 'T')] = {'A': configurator.lipids.A, 'B': 0.0, 'r0': configurator.lipids.r0, 'rc': configurator.lipids.rc} # Interface to tail is not sticky
        
        #Tail to other interactions
        glf.params[('T', 'T')] = {'A': configurator.lipids.A, 'B': configurator.lipids.B, 'r0': configurator.lipids.r0, 'rc': configurator.lipids.rc} # Tail to tail is sticky
        
        # Set the cutoff distance accordingly
        glf.r_cut[('H', 'H')] = configurator.lipids.rhh
        glf.r_cut[('H', 'I')] = configurator.lipids.r0
        glf.r_cut[('H', 'T')] = configurator.lipids.r0
        glf.r_cut[('I', 'I')] = configurator.lipids.rc
        glf.r_cut[('I', 'T')] = configurator.lipids.r0
        glf.r_cut[('T', 'T')] = configurator.lipids.rc
  
    # Only create AH domains if necessary
    if configurator.ahdomain:
        ###############################
        # AH-AH interactions
        ###############################
        # We also have to assign interactions between the ah pieces and every other body in the
        # system.  Self interactions for the ahs is not always zero!
        glf.params[('AH1', 'AH1')] = {'A': configurator.ahdomain.A, 'B': configurator.ahdomain.Bself_other, 'r0': configurator.ahdomain.r0, 'rc': configurator.ahdomain.rc}
        glf.params[('AH1', 'AH2')] = {'A': configurator.ahdomain.A, 'B': 0.0, 'r0': configurator.ahdomain.r0, 'rc': configurator.ahdomain.rc}
        glf.params[('AH2', 'AH2')] = {'A': configurator.ahdomain.A, 'B': configurator.ahdomain.Bself_hydrophobic, 'r0': configurator.ahdomain.r0, 'rc': configurator.ahdomain.rc}
        glf.r_cut[('AH1', 'AH1')] = configurator.ahdomain.rc
        glf.r_cut[('AH1', 'AH2')] = configurator.ahdomain.r0
        glf.r_cut[('AH2', 'AH2')] = configurator.ahdomain.rc
    
    # Create some simple sphere parameters
    if configurator.simple_spheres:
        ###############################
        # Simple sphere self-interactions
        ###############################
        glf.params[('SS1', 'SS1')] = {'A': configurator.simple_spheres.A, 'B': 0.0, 'r0': configurator.simple_spheres.r0, 'rc': configurator.simple_spheres.rc}
        glf.r_cut[('SS1', 'SS1')] = configurator.simple_spheres.rc
   
    if configurator.lipids and configurator.ahdomain:
        ###############################
        # Membrane-AH interactions
        ###############################
        # Initialize interactions to be purely repulsive to start with, use the average distance
        # as we have particles of different size interacting
        for x in itertools.product(lipid_types, ah_types):
            glf.params[x] = {'A': configurator.ahdomain.A, 'B': 0.0, 'r0': 0.5*(configurator.ahdomain.r0+configurator.lipids.r0), 'rc': 1.0*(configurator.ahdomain.r0+configurator.lipids.r0)}
            glf.r_cut[x] = 0.5*(configurator.ahdomain.r0+configurator.lipids.r0)
        
        # Set membrane-head ah-head interaction
        glf.params[('H', 'AH1')] = {'A': configurator.ahdomain.A, 'B': configurator.ahdomain.Bsurface, 'r0': 0.5*(configurator.ahdomain.r0+configurator.lipids.r0), 'rc': 1.0*(configurator.ahdomain.r0+configurator.lipids.r0)}
        glf.r_cut[('H', 'AH1')] = 1.0*(configurator.ahdomain.r0+configurator.lipids.r0)
        # Set the intermediate range interaction for AH-second
        glf.params[('I', 'AH2')] = {'A': configurator.ahdomain.A, 'B': configurator.ahdomain.Bintermediate, 'r0': 0.5*(configurator.ahdomain.r0+configurator.lipids.r0), 'rc': 1.0*(configurator.ahdomain.r0+configurator.lipids.r0)}
        glf.r_cut[('I', 'AH2')] = 1.0*(configurator.ahdomain.r0+configurator.lipids.r0)
        # Now set only the tail-tail interactions for ahs and the membrane
        glf.params[('T', 'AH2')] = {'A': configurator.ahdomain.A, 'B': configurator.ahdomain.Bdeep, 'r0': 0.5*(configurator.ahdomain.r0+configurator.lipids.r0), 'rc': 1.0*(configurator.ahdomain.r0+configurator.lipids.r0)}
        glf.r_cut[('T', 'AH2')] = 1.0*(configurator.ahdomain.r0+configurator.lipids.r0)

    if configurator.lipids and configurator.simple_spheres:
        ###############################
        # Membrane-Simple spheres interactions
        ###############################
        for x in itertools.product(lipid_types, ss_types):
            glf.params[x] = {'A': configurator.simple_spheres.A, 'B': 0.0, 'r0': 0.5*(configurator.lipids.r0+configurator.simple_spheres.r0), 'rc': 1.0*(configurator.lipids.r0+configurator.simple_spheres.r0)}
            glf.r_cut[x] = 1.0*(configurator.lipids.r0+configurator.simple_spheres.r0)

    if configurator.ahdomain and configurator.simple_spheres:
        ###############################
        # AHdomain-Simple spheres interactions
        ###############################
        for x in itertools.product(ah_types, ss_types):
            glf.params[x] = {'A': configurator.simple_spheres.A, 'B': 0.0, 'r0': 0.5*(configurator.ahdomain.r0+configurator.simple_spheres.r0), 'rc': 1.0*(configurator.ahdomain.r0+configurator.simple_spheres.r0)}
            glf.r_cut[x] = 1.0*(configurator.ahdomain.r0+configurator.simple_spheres.r0)
    
    ###############################
    # Bonded and angle interactions
    ###############################
    # Assign bonded interaction strengths
    if configurator.lipids or configurator.ahdomain:
        harmonic = md.bond.Harmonic()
        if configurator.lipids:
            harmonic.params['lipidbond'] = dict(k = configurator.lipids.kbond, r0 = configurator.lipids.rbond)
        if configurator.ahdomain:
            harmonic.params['ahbond'] = dict(k = configurator.ahdomain.kbond, r0 = configurator.ahdomain.rbond)
    
        # Assign angle interaction strengths
        angleharmonic = md.angle.Harmonic()
        if configurator.lipids:
            angleharmonic.params['lipidbend'] = dict(k = configurator.lipids.kbend, t0 = np.pi)
        if configurator.ahdomain:
            angleharmonic.params['ahbend'] = dict(k = configurator.ahdomain.kbend, t0 = configurator.ahdomain.thetabend)
    
    # Set up the integrator for the system
    integrator = md.Integrator(dt = configurator.deltatau)
    if configurator.integrator == 'langevin':
        langevin = hoomd.md.methods.Langevin(hoomd.filter.All(),
                                             kT = configurator.kT)
        if configurator.lipids:
            langevin.gamma['H'] = configurator.lipids.gamma_head
            langevin.gamma['I'] = configurator.lipids.gamma
            langevin.gamma['T'] = configurator.lipids.gamma
        if configurator.ahdomain:
            langevin.gamma['AH1'] = configurator.ahdomain.gamma
            langevin.gamma['AH2'] = configurator.ahdomain.gamma
        if configurator.simple_spheres:
            langevin.gamma['SS1'] = configurator.simple_spheres.gamma

        integrator.methods.append(langevin)

    elif configurator.integrator == 'NPT':
        npt = hoomd.md.methods.NPT(hoomd.filter.All(),
                                   kT = configurator.kT,
                                   tau = configurator.tau,
                                   S = 0.0,
                                   tauS = configurator.tauS,
                                   couple = "xy",
                                   box_dof = [True, True, False, False, False, False])

        integrator.methods.append(npt)

    elif configurator.integrator == "Brownian":
        brownian = hoomd.md.methods.Brownian(hoomd.filter.All(),
                                             kT = configurator.kT)
        if configurator.lipids:
            brownian.gamma['H'] = configurator.lipids.gamma_head
            brownian.gamma['I'] = configurator.lipids.gamma
            brownian.gamma['T'] = configurator.lipids.gamma
        if configurator.ahdomain:
            brownian.gamma['AH1'] = configurator.ahdomain.gamma
            brownian.gamma['AH2'] = configurator.ahdomain.gamma
        if configurator.simple_spheres:
            brownian.gamma['SS1'] = configurator.simple_spheres.gamma

        integrator.methods.append(brownian)
    
    
    integrator.forces.append(glf)
    if configurator.lipids or configurator.ahdomain:
        integrator.forces.append(harmonic)
        integrator.forces.append(angleharmonic)
    
    # Add to the simulation
    sim.operations.integrator = integrator
    
    # Run a simulation to get variables to work out, also thermalize the momenta
    sim.run(0)
    if configurator.init_type == 'all':
        if configurator.integrator == 'langevin':
            sim.state.thermalize_particle_momenta(hoomd.filter.All(), configurator.kT)
        elif configurator.integrator == 'NPT':
            sim.state.thermalize_particle_momenta(hoomd.filter.All(), configurator.kT)
            #npt.thermalize_thermostat_and_barostat_dof()
        elif configurator.integrator == 'Brownian':
            sim.state.thermalize_particle_momenta(hoomd.filter.All(), configurator.kT)
    elif configurator.init_type == 'read_gsd':
        if configurator.ahdomain:
            ahfilter = hoomd.filter.Type(['AH1', 'AH2'])
            sim.state.thermalize_particle_momenta(ahfilter, configurator.kT)
    
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
    
