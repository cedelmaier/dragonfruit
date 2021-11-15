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
    return N/A / (kc*q**4 + gamma*q**2)

# Radial average of points
def radial_average(uq, bin_size, pixel_size):
    cen_x = np.int32(uq.shape[1]/2)
    cen_y = np.int32(uq.shape[0]/2)
    X, Y = np.meshgrid(np.arange(uq.shape[1]) - cen_x, np.arange(uq.shape[0]) - cen_y)
    R = np.sqrt(X**2 + Y**2)
    Rreal = pixel_size * R

    # Rad keeps track of the index
    rad = np.arange(0, np.max(R), 1, dtype=np.int64)
    radreal = rad * pixel_size
    intensity = np.zeros(len(rad))
    #print(f"cen ({cen_x}, {cen_y})")
    #print(f"R: {R}")
    #print(f"Rreal: {Rreal}")
    #print(f"rad: {rad}")
    #print(f"radreal: {radreal}")
    index = 0
    for i in rad:
        currad = radreal[i]
        #print(f"{i} -> {currad}")
        mask = (np.greater_equal(Rreal, currad)) & (np.less(Rreal, currad + bin_size))
        #print(f"{mask}")
        values = np.abs(uq[mask])
        intensity[index] = np.mean(values)
        index += 1

    # Return the real radii that were evaluated and the corresponding 1D intensity
    return [radreal, intensity]


# Open the trajectory
#traj_all = gsd.hoomd.open('../data/traj_equil_membonly.gsd')
#traj_all = gsd.hoomd.open('../data/traj_postequil.gsd')
traj_all = gsd.hoomd.open('../data/20211102/strength_scan/traj_npt_membrane.gsd')
#traj_all = gsd.hoomd.open('../data/20211102/strength_scan/simulations/Bself1.166666_Bsurface2.33333_Bintermediate1.166666_Bdeep4.66666/s3/traj_langevin.gsd')
print(f"nframes: {len(traj_all)}")

# Set up some parameters for analysis
min_frame = 300
max_frame = 400
nmeasure = 10
xsize = 100.0
ysize = 100.0
Nx = 100
Ny = 100
Nxgrid = 100j
Nygrid = 100j
Ndirect = 20
deltaq = 0.0375
nlipids_per_layer = 10000

uq2d_fft = np.zeros((Nx, Ny, np.int32((max_frame - min_frame)/nmeasure)), dtype=np.complex128)
uq2d_direct = np.zeros((2*Ndirect+1, 2*Ndirect+1, np.int32((max_frame - min_frame)/nmeasure)), dtype=np.complex128)
area = np.zeros(np.int32((max_frame - min_frame)/nmeasure), dtype=np.float64)
qcutoff_arr = np.zeros(np.int32((max_frame - min_frame)/nmeasure), dtype=np.float64)

jdx = 0 # Measure which one we're doing
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
    Lx = traj.configuration.box[0]
    Ly = traj.configuration.box[1]
    area[jdx] = Lx*Ly

    # Figure out the membrane head indices
    [head_idx, leaf1_idx, leaf2_idx] = GetHeadIndices(traj, 2*nlipids_per_layer, 4)

    # Get the Z positions and re-center
    positions = traj.particles.position
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
    print(f"  sigma per pixel (FFT): {xpixel}, ypixel: {ypixel}") # sigma/pixel
    from scipy.interpolate import griddata
    grid_z1 = griddata(r1, z1, (grid_x, grid_y), method='cubic', fill_value = np.mean(z1))
    grid_z2 = griddata(r2, z2, (grid_x, grid_y), method='cubic', fill_value = np.mean(z2))

    # XXX: Maybe do something fancy with replacing the NaNs at the boundaries with a fill
    # value taken from a nearest neighbor type interpolation?

    ## Plot each if we so desire
    #plt.subplot(121)
    #plt.imshow(grid_z1.T, origin='lower')
    #plt.title('Z1')
    #plt.subplot(122)
    #plt.imshow(grid_z2.T, origin='lower')
    #plt.title('Z2')
    #plt.gcf().set_size_inches(6,6)
    #plt.show()

    # Take the 2d fourier transform
    u1 = np.fft.fft2(grid_z1, norm='forward')
    u1shift = np.fft.fftshift(u1)
    u2 = np.fft.fft2(grid_z2, norm='forward')
    u2shift = np.fft.fftshift(u2)

    # Get the sample frequencies
    freqx = np.fft.fftshift(np.fft.fftfreq(u1shift.shape[1], xpixel))
    freqy = np.fft.fftshift(np.fft.fftfreq(u1shift.shape[0], ypixel))
    fft_pixelsize = freqx[1] - freqx[0]

    # Sum the two now
    uq2d = 0.5*(u1shift + u2shift)
    uq2d_fft[:,:,jdx] = 0.5*(u1shift + u2shift)

    #plt.subplot(121)
    #plt.imshow(np.abs(uq2d.T), origin='lower')
    #plt.title('FFT')
    #plt.subplot(122)
    #plt.imshow(np.log10(np.square(np.abs(uq2d.T))), origin='lower')
    #plt.title('log10FFT squared')
    #plt.show()

    # Create the XY grid of the uq points for radial average
    [radii_fft, intensity_fft] = radial_average(uq2d, deltaq, fft_pixelsize)

    ################
    # Direct calculation with own wave-vector
    # Loop over wave vectors 2*pi*n/L
    u1direct = np.zeros((2*Ndirect+1,2*Ndirect+1), dtype=np.complex128)
    u2direct = np.zeros((2*Ndirect+1,2*Ndirect+1), dtype=np.complex128)

    for n in np.arange(-Ndirect,Ndirect+1,1):
        for m in np.arange(-Ndirect,Ndirect+1,1):
            idx = n + Ndirect
            kdx = m + Ndirect

            q = 2.0*np.pi*np.array([n/Lx, m/Ly])
            qnorm = np.linalg.norm(q)

            #print(f"    ({n}, {m})({idx}, {jdx}) -> ({q[0]}, {q[1]}) -> qnorm: {qnorm}")

            for k,r1k in enumerate(r1):
                val = z1[k] * np.exp(-1j*np.dot(q,r1k))
                u1direct[idx,kdx] += val
            for k,r2k in enumerate(r2):
                val = z2[k] * np.exp(-1j*np.dot(q,r2k))
                u2direct[idx,kdx] += val

    # Create an XY grid of the qvec points
    uq_direct = 1.0/(2.0*nlipids_per_layer)*(u1direct + u2direct)
    uq2d_direct[:,:,jdx] = uq_direct

    # Create the XY grid of the uq points for radial average
    [radii_direct, intensity_direct] = radial_average(uq_direct, deltaq, qcutoffx/2.0/np.pi)

    qcutoff_arr[jdx] = qcutoffx

    #plt.subplot(111)
    #plt.loglog(radii_fft, np.square(intensity_fft), 'r+')
    #plt.loglog(radii_direct, np.square(intensity_direct), 'bo')
    #plt.title('Step {} loglog'.format(itx))
    #plt.show()

    # Increment any counters we need
    jdx += 1

# What do we know about the entire simulation?
area_mean = np.mean(area)
qcutoff_mean = np.mean(qcutoff_arr)
print(f"Mean area (sigma): {area_mean}, mean APL: {area_mean/nlipids_per_layer}, qcutoff: {qcutoff_mean}")

# Now, compile all the results together into one plot
uq2d_fft_mean = np.mean(uq2d_fft, axis=2)
uq2d_direct_mean = np.mean(uq2d_direct, axis=2)
[radii_fft_mean, intensity_fft_mean] = radial_average(uq2d_fft_mean, deltaq, qcutoff_mean/2.0/np.pi)
[radii_direct_mean, intensity_direct_mean] = radial_average(uq2d_direct_mean, deltaq, qcutoff_mean/2.0/np.pi)

# Plot the data without fits
su_fft = np.square(intensity_fft_mean)
su_direct = np.square(intensity_direct_mean)
x_fft = radii_fft_mean
x_direct = radii_direct_mean

plt.style.use(plotstyle_stl)
fig, ax = plt.subplots(1, 1)
ax.scatter(x_fft, su_fft, color='r', marker='+', linewidth=1)
ax.scatter(x_direct, su_direct, color='b', marker='o', s=80, facecolors='none')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'q ($\sigma^{-1}$)')
ax.set_ylabel(r'$ \langle | u(q)^{2} | \rangle $ ($\sigma^{2}$)')
fig.tight_layout()
fig.savefig('suq_spectrum.pdf', dpi=fig.dpi)
plt.close(fig)

# Try to fit the data
# Fit the data using a lambda function
npoints_fit = 2
from scipy.optimize import curve_fit
#popt_fft, pcov_fft = curve_fit(lambda q, kc, gamma: suq_curve(q, nlipids_per_layer, area_mean, kc, gamma), x_fft[0:npoints_fit], su_fft[0:npoints_fit])
#print(f"Simulation fit values:")
#print(popt_fft)
#print(pcov_fft)
#popt_direct, pcov_direct = curve_fit(lambda q, kc, gamma: suq_curve(q, nlipids_per_layer, area_mean, kc, gamma), x_direct[0:npoints_fit], su_direct[0:npoints_fit])
#print(popt_direct)
#print(pcov_direct)
popt_fft, pcov_fft = curve_fit(lambda q, kc: suq_curve(q, nlipids_per_layer, area_mean, kc, 0.0), x_fft[0:npoints_fit], su_fft[0:npoints_fit])
print(f"Simulation fit values:")
print(popt_fft)
print(pcov_fft)
popt_direct, pcov_direct = curve_fit(lambda q, kc: suq_curve(q, nlipids_per_layer, area_mean, kc, 0.0), x_direct[0:npoints_fit], su_direct[0:npoints_fit])
print(popt_direct)
print(pcov_direct)

# Draw and save off the information
fig, ax = plt.subplots(1, 1)
#ax.plot(x_fft, suq_curve(x_fft, N = nlipids_per_layer, A = area_mean, kc= popt_fft[0], gamma = popt_fft[1]), 'k--')
#ax.plot(x_direct, suq_curve(x_direct, N = nlipids_per_layer, A = area_mean, kc = popt_direct[0], gamma = popt_direct[1]), 'k..')
ax.plot(x_fft, suq_curve(x_fft, N = nlipids_per_layer, A = area_mean, kc= popt_fft[0], gamma = 0.0), 'k--')
ax.plot(x_direct, suq_curve(x_direct, N = nlipids_per_layer, A = area_mean, kc = popt_direct[0], gamma = 0.0), 'k..')
ax.scatter(x_fft, su_fft, color='r', marker='+', linewidth=1)
ax.scatter(x_direct, su_direct, color='b', marker='o', s=80, facecolors='none')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'q ($\sigma^{-1}$)')
ax.set_ylabel(r'$ \langle | u(q)^{2} | \rangle $ ($\sigma^{2}$)')


fig.tight_layout()
fig.savefig('suq_spectrum_fit.pdf', dpi=fig.dpi)
plt.close(fig)

