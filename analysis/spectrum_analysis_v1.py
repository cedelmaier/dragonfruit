#!/usr/bin/env python3

import fresnel
import freud
import gsd.hoomd
import sys

import numpy as np
import matplotlib.pyplot as plt

# Set up a style for this
plotstyle_stl = {
    "axes.titlesize" : 18,
    "axes.labelsize": 16,
    "lines.linewidth" : 3,
    "lines.markersize" : 10,
    "xtick.labelsize" : 15,
    "ytick.labelsize" : 15,
    "font.size" : 15,
    "font.serif" : "Arial",
    "font.sans-serif" : "Arial"
    }

def GetHeadIndices(traj, nlipids, nbeads):
    nmembrane = nlipids * nbeads
    head_idx = np.arange(0, nmembrane, nbeads)
    idx_skip = nbeads*2
    leaf1_idx = np.arange(0, nmembrane, idx_skip)
    leaf2_idx = np.arange(nbeads, nmembrane, idx_skip)

    return [head_idx, leaf1_idx, leaf2_idx]

def suq_curve(q, N, A, kc, gamma):
    return (N/A) / (kc*q**4 + gamma*q**2)

# Radial average of points
def radial_average(uq, bin_size, pixel_size):
    cen_x = np.int32(uq.shape[1]/2)
    cen_y = np.int32(uq.shape[0]/2)
    X, Y = np.meshgrid(np.arange(uq.shape[1]) - cen_x, np.arange(uq.shape[0]) - cen_y)
    R = np.sqrt(X**2 + Y**2)
    # Rreal has the information about the actual location of the points
    Rreal = pixel_size * R
    #print(f"Rreal:\n{Rreal}")

    # Set an Rmax so that we don't go outside the circle averaging
    Rmax = Rreal[0, cen_x]
    #print(f"Rmax:\n{Rmax}")

    # Calculate the number of bins
    nbins = np.int32(np.floor(Rmax / bin_size))
    bin_edges = np.arange(0, nbins+1)*bin_size
    bin_mids = (bin_edges[0:-1] + bin_edges[1:])/2
    #print(f"nbins: {nbins}, bin_edges: {bin_edges}, bin_mids: {bin_mids}")

    intensity = np.zeros(len(bin_mids))
    # Loop over the bins and count up what we need to
    for ibin in np.arange(nbins):
        mask = (np.greater_equal(Rreal, bin_edges[ibin])) & (np.less(Rreal, bin_edges[ibin+1]))
        #print(f"{mask}")
        values = np.abs(uq[mask])
        #print(f"{values}")
        intensity[ibin] = np.mean(values)

    return [bin_mids, intensity]



# Open the trajectory
#traj_all = gsd.hoomd.open('../data/traj_equil_membonly.gsd')
#traj_all = gsd.hoomd.open('../data/traj_postequil.gsd')
traj_all = gsd.hoomd.open('./traj_npt_membrane.gsd')
#traj_all = gsd.hoomd.open('./traj_langevin_membrane.gsd')
#traj_all = gsd.hoomd.open('./traj_brownian_membrane.gsd')
#traj_all = gsd.hoomd.open('../data/20211102/strength_scan/simulations/Bself1.166666_Bsurface2.33333_Bintermediate1.166666_Bdeep4.66666/s3/traj_langevin.gsd')
print(f"nframes: {len(traj_all)}")

# Constants of conversion in our system
kB = 1.987204259e-3 # kcal/(mol * K)
kTroom = 0.5961 # 1 kT = 0.5961 kcal/mol
sigma = 0.75 # 0.75 nm per bead size, convert everything

# Set up some parameters for analysis
min_frame = 450
max_frame = 500
nmeasure = 10
Nx = 200
Ny = 200
Nxgrid = 200j
Nygrid = 200j
Ndirect = 3
deltaq = 0.05 # nm
nlipids_per_layer = 40000

uq_2d_fft        = np.zeros((Nx, Ny, np.int32((max_frame - min_frame)/nmeasure)), dtype=np.complex128)
uq_2d_direct    = np.zeros((2*Ndirect+1, 2*Ndirect+1, np.int32((max_frame - min_frame)/nmeasure)), dtype=np.complex128)
q_direct        = np.zeros((2*Ndirect+1, 2*Ndirect+1, np.int32((max_frame - min_frame)/nmeasure)), dtype=np.float64)
qcutoff_arr     = np.zeros(np.int32((max_frame - min_frame)/nmeasure), dtype=np.float64)
pixel_size_fft  = np.zeros(np.int32((max_frame - min_frame)/nmeasure), dtype=np.float64)
area            = np.zeros(np.int32((max_frame - min_frame)/nmeasure), dtype=np.float64)

jdx = 0 # Measure which one we're doing
Lmax = 0.0
# Loop over trajectories
for itx,traj in enumerate(traj_all):
    # Take some frames near the middle of the run
    if itx % nmeasure != 0:
        continue
    if itx < min_frame:
        continue
    if itx > max_frame:
        break

    print(f"Frame: {itx}")
    print(f"  nparticles: {traj.particles.N}")
    Lx = traj.configuration.box[0]*sigma
    Ly = traj.configuration.box[1]*sigma
    area[jdx] = Lx*Ly
    Lmax = np.maximum(Ly, np.maximum(Lx, Lmax))

    # Figure out the membrane head indices
    [head_idx, leaf1_idx, leaf2_idx] = GetHeadIndices(traj, 2*nlipids_per_layer, 4)

    # Get the Z positions and re-center
    positions = traj.particles.position*sigma
    z1 = positions[leaf1_idx, 2]
    z2 = positions[leaf2_idx, 2]
    z0 = (np.sum(z1) + np.sum(z2))/(len(z1) + len(z2))
    z1 = z1 - z0
    z2 = z2 - z0

    # Get the r positions for each leaflet
    r1 = positions[leaf1_idx, 0:2]
    r2 = positions[leaf2_idx, 0:2]

    # Other shared information
    qcutoffx = 2.0*np.pi/Lx
    qcutoffy = 2.0*np.pi/Ly
    print(f"  Lx: {Lx}, Ly: {Ly}")
    print(f"  Area: {area[jdx]}, APL: {area[jdx]/nlipids_per_layer}")
    print(f"  qcutoffx: {qcutoffx}, qcutoffy: {qcutoffy}")

    ################
    # Interpolation method to put on a grid
    # Create the grid for interpolation, taking care to remove NaN calculations
    grid_x, grid_y = np.mgrid[-Lx/2:Lx/2:Nxgrid,-Ly/2:Ly/2:Nygrid]
    # For the grid, calculate sigma/pixel
    xpixel = Lx/Nx
    ypixel = Ly/Ny
    print(f"  FFT:")
    print(f"    sample spacing (FFT): {xpixel}") # sigma/pixel
    from scipy.interpolate import griddata
    grid_z1 = griddata(r1, z1, (grid_x, grid_y), method='cubic', fill_value = np.mean(z1))
    grid_z2 = griddata(r2, z2, (grid_x, grid_y), method='cubic', fill_value = np.mean(z2))

    # Take the 2d fourier transform
    u1 = np.fft.fft2(grid_z1, norm='forward')
    u1shift = np.fft.fftshift(u1)
    u2 = np.fft.fft2(grid_z2, norm='forward')
    u2shift = np.fft.fftshift(u2)

    # Get the sample frequencies
    freqx = np.fft.fftshift(np.fft.fftfreq(u1shift.shape[1], xpixel))
    freqy = np.fft.fftshift(np.fft.fftfreq(u1shift.shape[0], ypixel))
    #fft_pixelsize = (freqx[1] - freqx[0])*(2.0*np.pi)
    fft_pixelsize = (freqx[1] - freqx[0])
    pixel_size_fft[jdx] = fft_pixelsize
    print(f"    frequency pixel size (FFT): {pixel_size_fft[jdx]}")

    # Sum the two now
    uq_2d_fft[:,:,jdx] = 0.5*(u1shift + u2shift)

    ################
    # Direct calculation with own wave-vector
    # Loop over wave vectors 2*pi*n/L
    print(f"  Direct:")
    u1direct = np.zeros((2*Ndirect+1,2*Ndirect+1), dtype=np.complex128)
    u2direct = np.zeros((2*Ndirect+1,2*Ndirect+1), dtype=np.complex128)
    qdirect  = np.zeros((2*Ndirect+1,2*Ndirect+1), dtype=np.float64)

    for n in np.arange(-Ndirect,Ndirect+1,1):
        for m in np.arange(-Ndirect,Ndirect+1,1):
            idx = n + Ndirect
            kdx = m + Ndirect

            q = 2.0*np.pi*np.array([n/Lx, m/Ly])
            qnorm = np.linalg.norm(q)

            print(f"    ({n}, {m})({idx}, {kdx}) -> ({q[0]}, {q[1]}) -> qnorm: {qnorm}")
            qdirect[idx,kdx] = qnorm

            for k,r1k in enumerate(r1):
                val = z1[k] * np.exp(-1j*np.dot(q,r1k))
                u1direct[idx,kdx] += val
            for k,r2k in enumerate(r2):
                val = z2[k] * np.exp(-1j*np.dot(q,r2k))
                u2direct[idx,kdx] += val

    uq_direct = 1.0/(2.0*nlipids_per_layer)*(u1direct + u2direct)

    # Create an XY grid of the qvec points
    uq_2d_direct[:,:,jdx] = uq_direct
    q_direct[:,:,jdx] = qdirect # Sampling frequency locations for direct measurements
    qcutoff_arr[jdx] = qcutoffx

    # Increment any counters we need
    jdx += 1

# What do we know about the entire simulation?
area_mean = np.mean(area)
qcutoff_mean = np.mean(qcutoff_arr)
print(f"Mean area (sigma): {area_mean}, mean APL: {area_mean/nlipids_per_layer}, qcutoff: {qcutoff_mean}")

## Do the averaging here, and pay attention to how it is done!
## XXX
## For each time point, break down the averaging
intensity_direct_list   = []
intensity_fft_list      = []
for itx in np.arange(uq_2d_direct.shape[-1]):
    print(f"Timeframe: {itx}")
    print(f"Direct")
    [radii_direct, intensity_direct] = radial_average(uq_2d_direct[:,:,itx], deltaq, qcutoff_arr[itx])
    #[radii_direct, intensity_direct] = radial_average(uq_2d_direct[:,:,itx], deltaq, qcutoff_arr[itx]/(2.0*np.pi))
    print(f"  radii_direct:       {radii_direct}")
    print(f"  intensity_direct:   {intensity_direct}")
    intensity_direct_list.append(intensity_direct)

    # Now do the FFT version
    print(f"FFT")
    print(f"  pixel_size: {pixel_size_fft[itx]}")
    [radii_fft, intensity_fft] = radial_average(uq_2d_fft[:,:,itx], deltaq, pixel_size_fft[itx]*(2.0*np.pi))
    #[radii_fft, intensity_fft] = radial_average(uq_2d_fft[:,:,itx], deltaq, pixel_size_fft[itx])
    print(f"  radii_fft:            {radii_fft}")
    print(f"  intensity_fft:        {intensity_fft}")
    intensity_fft_list.append(intensity_fft)

print(f"Averages:")
print(f"Direct list:    {intensity_direct_list}")
print(f"FFT list:       {intensity_fft_list}")
uq_direct_arr   = np.array(intensity_direct_list)
uq_fft_arr      = np.array(intensity_fft_list)
uq_direct_mean  = np.mean(uq_direct_arr, axis = 0)
uq_fft_mean     = np.mean(uq_fft_arr, axis = 0)
print(f"Direct mean:    {uq_direct_mean}")
print(f"FFT mean:       {uq_fft_mean}")
print(f"Cutoff frequency: {2.0*np.pi/Lmax}")

# Plot the data without fits
# NOTE: There is a huge thing here, as you need the number of lipids in the layer to scale correctly!!!!!
su_fft = np.square(uq_fft_mean)*nlipids_per_layer
su_direct = np.square(uq_direct_mean)*nlipids_per_layer
x_fft = radii_fft
x_direct = radii_direct

plt.style.use(plotstyle_stl)
fig, ax = plt.subplots(1, 1)
ax.scatter(x_fft, su_fft, color='r', marker='+', linewidth=1)
ax.scatter(x_direct, su_direct, color='b', marker='o', s=80, facecolors='none')
ax.axvline(x = qcutoff_mean, ymin = 0, ymax = 1.0, color = 'm', linestyle = '--')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'q ($nm^{-1}$)')
ax.set_ylabel(r'$ \langle | u(q) |^{2} \rangle $ ($nm^{2}$)')
fig.tight_layout()
fig.savefig('suq_spectrum.pdf', dpi=fig.dpi)
plt.close(fig)

# Try to fit the data
# Fit the data using a lambda function
npoints_fit = 10
from scipy.optimize import curve_fit
#popt_fft, pcov_fft = curve_fit(lambda q, kc, gamma: suq_curve(q, nlipids_per_layer, area_mean, kc, gamma), x_fft[0:npoints_fit], su_fft[0:npoints_fit])
#print(f"Simulation fit values:")
#print(popt_fft)
#print(pcov_fft)
#popt_direct, pcov_direct = curve_fit(lambda q, kc, gamma: suq_curve(q, nlipids_per_layer, area_mean, kc, gamma), x_direct[0:npoints_fit], su_direct[0:npoints_fit])
#print(popt_direct)
#print(pcov_direct)

# Generate a guess based on the first 2 points
kc_guess1 = nlipids_per_layer / area_mean / su_fft[1] / (x_fft[1]**4)
kc_guess2 = nlipids_per_layer / area_mean / su_fft[2] / (x_fft[2]**4)
print(f"x1: {x_fft[1]}, {su_fft[1]} -> kc: {kc_guess1}")
print(f"x2: {x_fft[2]}, {su_fft[2]} -> kc: {kc_guess2}")

#from scipy.optimize import fsolve
#kc_guess = fsolve(lambda kc: suq_curve(x_fft[0], nlipids_per_layer, area_mean, kc, 0.0) - su_fft[0], 0.01)
#print(f"kc guess: {kc_guess}")


popt_fft, pcov_fft = curve_fit(lambda q, kc: suq_curve(q, nlipids_per_layer, area_mean, kc, 0.0), x_fft[1:npoints_fit], su_fft[1:npoints_fit], kc_guess1)
print(f"Simulation fit values:")
print(popt_fft)
print(pcov_fft)
popt_direct, pcov_direct = curve_fit(lambda q, kc: suq_curve(q, nlipids_per_layer, area_mean, kc, 0.0), x_direct[1:npoints_fit], su_direct[1:npoints_fit], kc_guess2)
print(popt_direct)
print(pcov_direct)

# Draw and save off the information
fig, ax = plt.subplots(1, 1)
#ax.plot(x_fft, suq_curve(x_fft, N = nlipids_per_layer, A = area_mean, kc= popt_fft[0], gamma = popt_fft[1]), 'k--')
#ax.plot(x_direct, suq_curve(x_direct, N = nlipids_per_layer, A = area_mean, kc = popt_direct[0], gamma = popt_direct[1]), 'k..')
ax.plot(x_fft[1:], suq_curve(x_fft[1:], N = nlipids_per_layer, A = area_mean, kc = popt_fft[0], gamma = 0.0), 'k--')
ax.plot(x_direct[1:], suq_curve(x_direct[1:], N = nlipids_per_layer, A = area_mean, kc = popt_direct[0], gamma = 0.0), 'k:')
#ax.plot(x_fft[1:], suq_curve(x_fft[1:], N = nlipids_per_layer, A = area_mean, kc = 150/4.1, gamma = 0.0), 'm-.')
ax.scatter(x_fft, su_fft, color='r', marker='+', linewidth=1)
ax.scatter(x_direct, su_direct, color='b', marker='o', s=80, facecolors='none')
ax.axvline(x = qcutoff_mean, ymin = 0, ymax = 1.0, color = 'm', linestyle = '--')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'q ($nm^{-1}$)')
ax.set_ylabel(r'$ \langle | u(q) |^{2} \rangle $ ($nm^{2}$)')


fig.tight_layout()
fig.savefig('suq_spectrum_fit.pdf', dpi=fig.dpi)
plt.close(fig)

