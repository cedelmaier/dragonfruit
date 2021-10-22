#!/usr/bin/env python3

import fresnel
import freud
import gsd.hoomd
import sys

import numpy as np
import matplotlib.pyplot as plt

def GetHeadIndices(traj, nlipids, nbeads):
    nmembrane = nlipids * nbeads
    head_idx = np.arange(0, nmembrane, nbeads)
    idx_skip = nbeads*2
    leaf1_idx = np.arange(0, nmembrane, idx_skip)
    leaf2_idx = np.arange(nbeads, nmembrane, idx_skip)

    return [head_idx, leaf1_idx, leaf2_idx]

# Open the trajectory
traj_all = gsd.hoomd.open('../data/traj_equil_membonly.gsd')
#traj_all = gsd.hoomd.open('../data/traj_postequil.gsd')

# Loop over trajectories
uq_abs2mean = np.zeros(72)
suq_mean = np.zeros(101)
ntraj = len(traj_all)
print(ntraj)
for itx,traj in enumerate(traj_all):
    # Take some frames near the middle of the run
    if itx % 10 != 0:
        continue
    if itx < 500:
        continue
    if itx >= 600:
        break

    print(itx)
    Lx = traj.configuration.box[0]
    Ly = traj.configuration.box[1]

    # Figure out the membrane head indices
    nlipids_per_layer = 10000
    [head_idx, leaf1_idx, leaf2_idx] = GetHeadIndices(traj, 20000, 4)

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
    print(f"Lx: {Lx}, Ly: {Ly}")
    print(f"qcutoffx: {qcutoffx}, qcutoffy: {qcutoffy}")

    #################
    ## Interpolation method to put on a grid
    ## Create the grid for interpolation, taking care to remove NaN calculations
    #xsize = 100
    #ysize = 100
    #grid_x, grid_y = np.mgrid[-xsize/2:xsize/2:101j,-ysize/2:ysize/2:101j]
    #from scipy.interpolate import griddata
    #grid_z1 = griddata(r1, z1, (grid_x, grid_y), method='cubic')
    #grid_z2 = griddata(r2, z2, (grid_x, grid_y), method='cubic')

    ## Remove NaN calculations from interpolation
    #grid_z1[np.isnan(grid_z1)] = np.mean(z1)
    #grid_z2[np.isnan(grid_z2)] = np.mean(z2)

    ## Plot each if we so desire
    #plt.imshow(grid_z1.T, extent=(-Lx/2,Lx/2,-Ly/2,Ly/2), origin='lower')
    #plt.show()
    #plt.imshow(grid_z2.T, extent=(-Lx/2,Lx/2,-Ly/2,Ly/2), origin='lower')
    #plt.show()

    ## Take the 2d fourier transform
    #u1 = np.fft.fft2(grid_z1)
    #u2 = np.fft.fft2(grid_z2)

    ## Compute the 2d combined value
    #uq2d = 1.0/(2.0*nlipids_per_layer)*(u1 + u2)
    #plt.imshow(np.abs(uq2d.T), origin='lower')
    #plt.show()
    #plt.imshow(np.log10(np.square(np.abs(uq2d.T))), origin='lower')
    #plt.show()

    ## Create the XY grid of the uq points for radial average
    #cen_x = uq2d.shape[1]/2+1
    #cen_y = uq2d.shape[0]/2+1
    #X,Y = np.meshgrid(np.arange(uq2d.shape[1]) - cen_x, np.arange(uq2d.shape[0]) - cen_y)
    #R = np.sqrt(X**2+Y**2)
    #rad = np.arange(1, np.max(R), 1)
    #intensity = np.zeros(len(rad))
    #index = 0
    #bin_size = 1
    #for i in rad:
    #    mask = (np.greater(R, i - bin_size) & np.less(R, i + bin_size))
    #    values = np.abs(uq2d[mask])
    #    intensity[index] = np.mean(values)
    #    index += 1

    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(rad, intensity, linewidth=2)
    #plt.show()

    ## Add to the main uq_mean
    #uq_abs2mean += np.square(intensity)

    #sys.exit(1)

    ################
    # Direct calculation with own wave-vector
    # Loop over wave vectors 2*pi*n/L
    Nx = 10
    #u1direct = np.zeros((2*Nx+1,2*Nx+1), dtype=np.complex128)
    #u2direct = np.zeros((2*Nx+1,2*Nx+1), dtype=np.complex128)

    # Try directly creating the linear versions as we go and bin them appropriately
    deltaq = 0.05
    qvar = np.linspace(0, 100*deltaq, 101)
    u1linear = np.zeros(101, dtype=np.complex128)
    u1nbin = np.zeros(101, dtype=np.int32)
    u2linear = np.zeros(101, dtype=np.complex128)
    u2nbin = np.zeros(101, dtype=np.int32)

    #qx, qy = np.meshgrid(2.0*np.pi/Lx*np.arange(-Nx,Nx+1,1), 2.0*np.pi/Ly*np.arange(-Nx,Nx+1,1))

    for n in np.arange(-Nx,Nx+1,1):
        for m in np.arange(-Nx,Nx+1,1):
            idx = n + Nx
            jdx = m + Nx

            q = 2.0*np.pi*np.array([n/Lx, m/Ly])
            qnorm = np.linalg.norm(q)
            # Figure out the linear bin this should be in
            linear_bin = np.int32(np.floor(qnorm/deltaq))

            print(f"({n}, {m}) -> ({q[0]}, {q[1]}) -> qnorm: {qnorm}, bin: {linear_bin}")

            for k,r1k in enumerate(r1):
                val = z1[k] * np.exp(-1j*np.dot(q,r1k))
                #u1direct[idx,jdx] += val
                u1linear[linear_bin] += val
                u1nbin[linear_bin] += 1
            for k,r2k in enumerate(r2):
                val = z2[k] * np.exp(-1j*np.dot(q,r2k))
                #u2direct[idx,jdx] += val
                u2linear[linear_bin] += val
                u2nbin[linear_bin] += 1

    # Create mean of the linear arrays
    u1lin = np.divide(u1linear, u1nbin)
    u1lin[np.isnan(u1lin)] = 0.0
    u2lin = np.divide(u2linear, u2nbin)
    u1lin[np.isnan(u2lin)] = 0.0
    ulin = 1.0/(2.0*nlipids_per_layer)*(u1lin+u2lin)
    suq = np.square(np.abs(ulin))
    suq_mean += suq

    #plt.scatter(x = qvar, y = suq, marker = 's', color = 'b')
    #plt.show()

    ## Create an XY grid of the qvec points
    #uq = 1.0/(2.0*nlipids_per_layer)*(u1direct + u2direct)

    #cen_x = uq.shape[1]/2+1
    #cen_y = uq.shape[0]/2+1
    #X,Y = np.meshgrid(np.arange(uq.shape[1]) - cen_x, np.arange(uq.shape[0]) - cen_y)
    #R = np.sqrt(X**2+Y**2)
    #rad = np.arange(1, np.max(R), 1)
    #intensity = np.zeros(len(rad))
    #index = 0
    ##bin_size = 0.7957747155 # Bin size in pixels to wave-vector
    #bin_size = 1
    #for i in rad:
    #    mask = (np.greater(R, i - bin_size) & np.less(R, i + bin_size))
    #    values = np.abs(uq[mask])
    #    intensity[index] = np.mean(values)
    #    index += 1

    #plt.plot(rad*2.0*np.pi/Lx, np.square(intensity))
    #plt.show()

# Calculate the final numbers
suq_ensemble = suq_mean/(itx+1)
plt.scatter(x = qvar, y = suq_ensemble, marker = 's', color = 'b')
plt.show()

## Calculate the average undulation spectrum
#suq = uq_abs2mean/itx
#    
#fig = plt.figure()
#plt.loglog(rad, suq)
#plt.show()
