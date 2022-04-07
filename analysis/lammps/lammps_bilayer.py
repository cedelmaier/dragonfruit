#!/usr/bin/env python3

import os
import re
import sys

import freud
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Magic to get other definitions in place
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'lib'))
from common import radial_average, ragged_mean
from stylelib.common_styles import septin_runs_stl

def Compute_U_FFT(Lx, Nx, r, z):
    xpixel = Lx/Nx
    Nxgrid = Nx * 1j

    grid_x, grid_y = np.mgrid[-Lx/2:Lx/2:Nxgrid, -Lx/2:Lx/2:Nxgrid]
    from scipy.interpolate import griddata
    grid_z = griddata(r, z, (grid_x, grid_y), method = 'cubic')
    grid_z_nearest = griddata(r, z, (grid_x, grid_y), method = 'nearest')
    grid_z[np.isnan(grid_z)] = grid_z_nearest[np.isnan(grid_z)]

    u = np.fft.fft2(grid_z) / (Nx * Nx)
    ushift = np.fft.fftshift(u)
    freqx = np.fft.fftshift(np.fft.fftfreq(ushift.shape[1], xpixel))

    uq_2d_fft_qcutoff = (freqx[1] - freqx[0])
    return [ushift, uq_2d_fft_qcutoff]

def Compute_U_DirectFast(Lx, ndirect, r, z):
    xvec = 2.0*np.pi/Lx*np.linspace(-ndirect, ndirect, 2*ndirect+1)
    yvec = 2.0*np.pi/Lx*np.linspace(-ndirect, ndirect, 2*ndirect+1)
    qmesh = np.array(np.meshgrid(xvec, yvec)).T.reshape(-1,2)
    udirect_fast = (np.sum(z * np.exp(-1j*np.dot(qmesh, r.T)), axis=-1)).reshape(2*ndirect+1,2*ndirect+1)

    return udirect_fast

def suq_curve(q, N, A, kc, gamma):
    return (N/A) / (kc*q**4 + gamma*q**2)

#datafile = "/Users/cedelmaier/Projects/Biophysics/septin_project/supra_cg/dragonfruit/data/20220311/lammps_hoomd_comparison/lammps/blm.lip.lammpstrj"
#datafile = "/Users/cedelmaier/Projects/Biophysics/septin_project/supra_cg/dragonfruit/data/20220321/lammps_hoomd_comparison/lammps/nph_langevin/blm.lip.lammpstrj"
#datafile = "/Users/cedelmaier/Projects/Biophysics/septin_project/supra_cg/dragonfruit/data/20220321/lammps_hoomd_comparison/lammps/nvt_langevin/blm.lip.lammpstrj"
#psffile  = "/Users/cedelmaier/Projects/Biophysics/septin_project/supra_cg/dragonfruit/data/20220311/lammps_hoomd_comparison/lammps/myfile.psf"

datafile = sys.argv[1]

# Unfortunately, we have to make our own file reader for the LAMMPS trajectory
natoms = -1
timestep = -1
box_data = np.zeros((3, 2), dtype = np.float32)

plt.style.use(septin_runs_stl)

# Actual freud data
box = None
data = None
modedata = {}
nframes = -1

min_frame = 404
#min_frame = 494
#min_frame = 502
max_frame = 504

timesteps = []
Nx = 200
Nxgrid = Nx * 1j
Ndirect = 1

with open(datafile, "r") as stream:
    for index, line in enumerate(stream):
        rline = line.rstrip()
        #print(rline)
        # Look for a timestep
        if rline == "ITEM: TIMESTEP":
            timestep = np.int32(next(stream).rstrip())

        # Look for a number of atoms
        if rline == "ITEM: NUMBER OF ATOMS":
            natoms = np.int32(next(stream).rstrip())

        # Look for the box bounds
        if rline == "ITEM: BOX BOUNDS pp pp pp":
            xlim = next(stream).rstrip().split()
            ylim = next(stream).rstrip().split()
            zlim = next(stream).rstrip().split()

            box_data[0,0] = np.float32(xlim[0])
            box_data[0,1] = np.float32(xlim[1])
            box_data[1,0] = np.float32(ylim[0])
            box_data[1,1] = np.float32(ylim[1])
            box_data[2,0] = np.float32(zlim[0])
            box_data[2,1] = np.float32(zlim[1])

            box = freud.box.Box.from_box(box_data[:, 1] - box_data[:, 0])

        # Look for the atoms
        if rline == "ITEM: ATOMS id type mol x y z":
            nframes += 1
            # Check for an early exit on the frame reading
            if nframes < min_frame:
                continue
            if nframes >= max_frame:
                break
            # Now we can read in the particles and their positions
            data = np.zeros((natoms,6), dtype = np.float32)
            timesteps.append(timestep)

            for idx in range(natoms):
                aline = next(stream).rstrip().split()

                data[idx,0] = np.float32(aline[0])
                data[idx,1] = np.float32(aline[1])
                data[idx,2] = np.float32(aline[2])
                data[idx,3] = np.float32(aline[3])
                data[idx,4] = np.float32(aline[4])
                data[idx,5] = np.float32(aline[5])

                #print(f"IDX: {idx}, type: {data[idx,1]}, position: {data[idx,3]}, {data[idx,4]}, {data[idx,5]}")
            nlipids_per_leaflet = natoms/4/2
            print(f"Timestep: {timestep}")
            print(f"Natoms: {natoms}")
            print(f"box: {box}")
            #print(f"data: {data}")
            print(f"nlipids_per_leaflet: {nlipids_per_leaflet}")
            print(f"frame: {nframes}")

            # Now that we have the data, process it! Do the calculation here, then just let it keep
            # looping on without changing anything.
            # Get the head indicies
            leaf1_h_idx = np.int32(np.where((data[:,1] == 1) & (data[:,0] < natoms/2))).flatten()
            leaf2_h_idx = np.int32(np.where((data[:,1] == 1) & (data[:,0] >= natoms/2))).flatten()
            #print(leaf1_h_idx)
            #print(leaf2_h_idx)

            Lx = box.Lx
            Ly = box.Ly
            qcutoffx = 2.0*np.pi/Lx

            positions = data[:,3:]
            z1 = positions[leaf1_h_idx, 2]
            z2 = positions[leaf2_h_idx, 2]
            z0 = (np.sum(z1) + np.sum(z2))/(len(z1) + len(z2))
            z1 = z1 - z0
            z2 = z2 - z0

            r1 = positions[leaf1_h_idx, 0:2]
            r2 = positions[leaf2_h_idx, 0:2]

            # Interpolation
            [ushift1, qcutoff1] = Compute_U_FFT(Lx, Nx, r1, z1)
            [ushift2, qcutoff2] = Compute_U_FFT(Lx, Nx, r2, z2)
            uq_2d_fft = 0.5*(ushift1 + ushift2)

            ## Direct fast measurement
            #udirectfast1 = Compute_U_DirectFast(Lx, Ndirect, r1, z1)
            #udirectfast2 = Compute_U_DirectFast(Lx, Ndirect, r2, z1)
            #uq_2d_direct_fast = 1.0/(2.0*nlipids_per_leaflet)*(udirectfast1 + udirectfast2)

            # Save off information for later!
            if 'uq_2d_fft_modes' not in modedata:
                modedata['uq_2d_fft_modes'] = {}
            if 'uq_2d_fft_qcutoff' not in modedata:
                modedata['uq_2d_fft_qcutoff'] = {}
            #if 'uq_2d_direct_modes' not in modedata:
            #    modedata['uq_2d_direct_modes'] = {}
            #if 'uq_2d_direct_qcutoff' not in modedata:
            #    modedata['uq_2d_direct_qcutoff'] = {}
            if 'area' not in modedata:
                modedata['area'] = {}

            modedata['uq_2d_fft_modes'][timestep] = uq_2d_fft
            modedata['uq_2d_fft_qcutoff'][timestep] = qcutoff1*2.0*np.pi

            #modedata['uq_2d_direct_modes'][timestep] = uq_2d_direct_fast
            #modedata['uq_2d_direct_qcutoff'][timestep] = qcutoffx

            modedata['area'][timestep] = Lx*Ly


# In theory, we now have all the membrane mode information, process it
print(f"Finished processing, now running analysis")
nframes_calculated = len(modedata['uq_2d_fft_modes'])
print(f"Number of frames: {nframes_calculated}")

deltaq = 0.0375

# Timepoints
timestep = np.array(timesteps, dtype=np.int64)
# Area
area_list               = modedata['area']
area_arr                = np.array([area_list[ts] for ts in timestep], dtype=np.float64)
# FFT calculation
uq_2d_fft_modes         = modedata['uq_2d_fft_modes']
uq_2d_fft_modes_arr     = np.array([uq_2d_fft_modes[ts] for ts in timestep], dtype=np.complex128)
uq_2d_fft_qcutoff       = modedata['uq_2d_fft_qcutoff']
uq_2d_fft_qcutoff_arr   = np.array([uq_2d_fft_qcutoff[ts] for ts in timestep], dtype=np.float64)
## Direct modes?
#uq_2d_direct_modes          = modedata['uq_2d_direct_modes']
#uq_2d_direct_modes_arr      = np.array([uq_2d_direct_modes[ts] for ts in timestep], dtype = np.complex128)
#uq_2d_direct_qcutoff        = modedata['uq_2d_direct_qcutoff']
#uq_2d_direct_qcutoff_arr    = np.array([uq_2d_direct_qcutoff[ts] for ts in timestep], dtype = np.float64)

# Loop over membrane modesl to calculate the max size of arrays
max_len = 0
for itx in np.arange(uq_2d_fft_modes_arr.shape[0]):
    [radii_fft, intensity_fft] = radial_average(uq_2d_fft_modes_arr[itx,:,:], deltaq, uq_2d_fft_qcutoff_arr[itx])
    if len(intensity_fft) > max_len:
        max_len = len(intensity_fft)
# Compute just like the real version
radii_fft_list = []
#radii_direct_list = []
intensity_fft_list = []
#intensity_direct_list = []
uq_2d_fft_qcutoff_list = []
#uq_2d_direct_qcutoff_list = []
area_list = []
for itx in np.arange(uq_2d_fft_modes_arr.shape[0]):
    [radii_fft, intensity_fft] = radial_average(uq_2d_fft_modes_arr[itx,:,:], deltaq, uq_2d_fft_qcutoff_arr[itx])
    intensity_fft_list.append(intensity_fft)
    radii_fft_list.append(radii_fft)
    uq_2d_fft_qcutoff_list.append(uq_2d_fft_qcutoff_arr[itx])
    
    #[radii_direct, intensity_direct] = radial_average(uq_2d_direct_modes_arr[itx,:,:], deltaq, uq_2d_direct_qcutoff_arr[itx])
    #intensity_direct_list.append(intensity_direct)
    #radii_direct_list.append(radii_direct)
    #uq_2d_direct_qcutoff_list.append(uq_2d_direct_qcutoff_arr[itx])

    area_list.append(area_arr[itx])

[radii_fft_mean, radii_fft_std] = ragged_mean(radii_fft_list)
[intensity_fft_mean, intensity_fft_std] = ragged_mean(intensity_fft_list)
su_fft = np.square(intensity_fft_mean)*nlipids_per_leaflet

#[radii_direct_mean, radii_direct_std] = ragged_mean(radii_direct_list)
#[intensity_direct_mean, intensity_direct_std] = ragged_mean(intensity_direct_list)
#su_direct = np.square(intensity_direct_mean)*nlipids_per_leaflet

area_mean = np.mean(area_list)

fig, ax = plt.subplots(1, 1, figsize=(15,10))

# Plot everything
ax.scatter(radii_fft_mean[1:], su_fft[1:], color = 'b', marker = '+', linewidth = 1)
#ax.scatter(radii_direct_mean[1:], su_direct[1:], color = 'r', marker = 'o', s = 80, facecolors = 'none')

# Figure out where cutoffs are
qcutoff_mean = np.mean(uq_2d_fft_qcutoff_list)
print(f"qcutoff_mean    = {qcutoff_mean}")
print(f"area_mean       = {area_mean}")
idx = np.where(np.greater(radii_fft_mean, qcutoff_mean))
idx = np.int32(idx[0][0])
jdx = np.where(np.greater(radii_fft_mean, 1.0))
jdx = np.int32(jdx[0][0])
# Add 1 to idx to correct pathological behavior?
#idx += 1

# Generate a guess
kcguess1 = 1.0*nlipids_per_leaflet / area_mean / su_fft[idx] / (radii_fft_mean[idx]**4)
# Try the fits
from scipy.optimize import curve_fit
popt_fft_kc, pcov_fft_kc = curve_fit(lambda q, kc: suq_curve(q, nlipids_per_leaflet, area_mean, kc, 0.0), radii_fft_mean[idx:jdx], su_fft[idx:jdx], bounds = ([0.0, np.inf]), p0 = [kcguess1])
popt_fft_ga, pcov_direct_gc = curve_fit(lambda q, kc, gamma: suq_curve(q, nlipids_per_leaflet, area_mean, kc, gamma), radii_fft_mean[idx:jdx], su_fft[idx:jdx], bounds = ([0.0, -np.inf], [np.inf, np.inf]), p0 = [kcguess1, 0.0])
#popt_direct_kc, pcov_direct_kcA = curve_fit(lambda q, kc: suq_curve(q, nlipids_per_leaflet, area_mean, kc, 0.0), radii_direct_mean[idx:jdx], su_direct[idx:jdx], bounds = ([0.0, np.inf]), p0 = [kcguess1])
#popt_direct_ga, pcov_direct_ga = curve_fit(lambda q, kc, gamma: suq_curve(q, nlipids_per_leaflet, area_mean, kc, gamma), radii_direct_mean[idx:jdx], su_direct[idx:jdx], bounds = ([0.0, -np.inf], [np.inf, np.inf]), p0 = [kcguess1, 0.0])

print(f"Simulation fit values:")
print(f"  kc(guess)         = {kcguess1}")
print(f"----No gamma----")
print(f"  FFT kc                = {popt_fft_kc[0]}")
#print(f"  Direct kc             = {popt_direct_kc[0]}")
print(f"---With gamma----")
print(f"  FFT kc, gamma         = {popt_fft_ga[0]}, {popt_fft_ga[1]}")
#print(f"  Direct kc, gamma      = {popt_direct_ga[0]}, {popt_direct_ga[1]}")

ax.plot(radii_fft_mean[idx:jdx], suq_curve(radii_fft_mean[idx:jdx], N = nlipids_per_leaflet, A = area_mean, kc = popt_fft_kc[0], gamma = 0.0), color = 'b', linestyle = '--')
ax.plot(radii_fft_mean[idx:jdx], suq_curve(radii_fft_mean[idx:jdx], N = nlipids_per_leaflet, A = area_mean, kc = popt_fft_ga[0], gamma = popt_fft_ga[1]), color = 'b', linestyle = ':')
#ax.plot(radii_direct_mean[idx:jdx], suq_curve(radii_direct_mean[idx:jdx], N = nlipids_per_leaflet, A = area_mean, kc = popt_direct_kc[0], gamma = 0.0), color = 'r', linestyle = '--')
#ax.plot(radii_direct_mean[idx:jdx], suq_curve(radii_direct_mean[idx:jdx], N = nlipids_per_leaflet, A = area_mean, kc = popt_direct_ga[0], gamma = popt_direct_ga[1]), color = 'r', linestyle = ':')

# Plot the cutoff
ax.axvline(x = qcutoff_mean, ymin = 0, ymax = 1.0, color = 'k', linestyle = '-')


# Set the log scale stuff
ax.set_ylim(1e-1,1e6)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_title('Membrane Modes')
ax.set_xlabel(r'q ($\sigma^{-1}$)')
ax.set_ylabel(r'$ N \langle | u(q) |^{2} \rangle $ ($\sigma^{2}$)')

fig.tight_layout()
fig.savefig('lammps_membranemodes.pdf', dpi = fig.dpi)

# Create a dataframe of the results so that we can easily dump to a CSV file
dfs = []
df_x_fft        = pd.DataFrame(radii_fft_mean, columns = ['x_fft'])
df_su_fft       = pd.DataFrame(su_fft, columns = ['su_fft'])
#df_x_direct     = pd.DataFrame(radii_direct_mean, columns = ['x_direct'])
#df_su_direct    = pd.DataFrame(su_direct, columns = ['su_direct'])
df_other        = pd.DataFrame([area_mean, nlipids_per_leaflet], columns = ['other'])
dfs.append(df_x_fft)
dfs.append(df_su_fft)
#dfs.append(df_x_direct)
#dfs.append(df_su_direct)
dfs.append(df_other)

# Combine all together
df = pd.concat(dfs, axis=1)

# Write to a dumpfile
with open('dumpdata.csv', 'w') as stream:
    df.to_csv(stream, index=False)

