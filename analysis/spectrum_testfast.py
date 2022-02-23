#!/usr/bin/env python3

import sys

import numpy as np
import matplotlib.pyplot as plt

import scipy
import timeit

# Defined function for use
def func(x, y):
    return np.cos(4*np.pi*x/100) * np.sin(8*np.pi*y/100)

# Define a radial average
def radial_average(uq, bin_size, pixel_size):
    cen_x = np.int32(uq.shape[1]/2)
    cen_y = np.int32(uq.shape[0]/2)
    X, Y = np.meshgrid(np.arange(uq.shape[1]) - cen_x, np.arange(uq.shape[0]) - cen_y)
    R = np.sqrt(X**2 + Y**2)
    Rreal = pixel_size * R

    Rmax = Rreal[0, cen_x]

    nbins = np.int32(np.floor(Rmax / bin_size))
    bin_edges = np.arange(0, nbins+1)*bin_size
    bin_mids = (bin_edges[0:-1] + bin_edges[1:])/2

    intensity = np.zeros(len(bin_mids))
    for ibin in np.arange(nbins):
        mask = (np.greater_equal(Rreal, bin_edges[ibin])) & (np.less(Rreal, bin_edges[ibin+1]))
        values = np.abs(uq[mask])
        intensity[ibin] = np.mean(values)

    return [bin_mids, intensity]

# Compute the FFT version
def compute_u_fft(Lx, Ly, Nx, Ny, r, z):
    xpixel = Lx/Nx
    ypixel = Ly/Ny
    Nxgrid = Nx * 1j
    Nygrid = Ny * 1j
    grid_x, grid_y = np.mgrid[-Lx/2:Lx/2:Nxgrid, -Ly/2:Ly/2:Nygrid]
    from scipy.interpolate import griddata
    grid_z = griddata(r, z, (grid_x, grid_y), method = 'cubic')
    grid_z_nearest = griddata(r, z, (grid_x, grid_y), method = 'nearest')
    grid_z[np.isnan(grid_z)] = grid_z_nearest[np.isnan(grid_z)]

    # Take the 2D fourier transform
    u = np.fft.fft2(grid_z, norm = 'forward')
    ushift = np.fft.fftshift(u)
    freqx = np.fft.fftshift(np.fft.fftfreq(ushift.shape[1], xpixel))
    freqy = np.fft.fftshift(np.fft.fftfreq(ushift.shape[0], ypixel))

    # Save the information
    u_fft = ushift
    qcutoff_fft = (freqx[1] - freqx[0])*2.0*np.pi

    return [u_fft, qcutoff_fft]

# Compute the slow version of the FFT
def compute_u_direct_slow(Lx, Ly, ndirect, r, z):
    udirect = np.zeros((2*ndirect+1,2*ndirect+1), dtype=np.complex128)
    
    for n in np.arange(-ndirect, ndirect+1, 1):
        for m in np.arange(-ndirect, ndirect+1, 1):
            idx = n + ndirect
            kdx = m + ndirect
    
            q = 2.0*np.pi*np.array([n/Lx, m/Ly])
    
            for k,rk in enumerate(r):
                val = z[k] * np.exp(-1j*np.dot(q,rk))
                udirect[idx,kdx] += val

    udirect = udirect / len(z)
    return [udirect]

# Fast calculation of direct
def compute_u_direct_fast(Lx, Ly, ndirect, r, z):
    # Create a mesh-based verson for fast array stuff
    xvec = 2.0*np.pi/Lx*np.linspace(-ndirect,ndirect,2*ndirect+1)
    yvec = 2.0*np.pi/Ly*np.linspace(-ndirect,ndirect,2*ndirect+1)
    qmesh = np.array(np.meshgrid(xvec, yvec)).T.reshape(-1,2)
    udirect_fast = (np.sum(z * np.exp(-1j*np.dot(qmesh, r.T)), axis=-1)).reshape(2*ndirect+1,2*ndirect+1) / len(z) 

    return [udirect_fast]

# Compute FFT timing
def fft_time():
    SETUP_CODE = '''
from __main__ import compute_u_fft, func
import numpy as np
'''
    
    TEST_CODE = '''
Nx = 200
Ny = 200
Lx = 200
Ly = 200
npoints = 2000
Ndirect = 21

rng = np.random.default_rng()
points = (rng.random((npoints, 2)) - 0.5)*(Lx)
values = func(points[:,0], points[:,1])
values = np.array([values])
positions = np.concatenate((points, values.T), axis = 1)
z = positions[:,2]
r = positions[:,0:2]
[u_fft, qcutoff_fft] = compute_u_fft(Lx, Ly, Nx, Ny, r, z)'''

    # timeit.repeat statement
    times = timeit.repeat(setup = SETUP_CODE,
                          stmt = TEST_CODE,
                          repeat = 3,
                          number = 100)
    print('FFT time:  {}'.format(min(times)))

# Compute direct(slow) time
def slow_time():
    SETUP_CODE = '''
from __main__ import compute_u_direct_slow, func
import numpy as np
'''

    TEST_CODE = '''
Nx = 200
Ny = 200
Lx = 200
Ly = 200
npoints = 2000
Ndirect = 5

rng = np.random.default_rng()
points = (rng.random((npoints, 2)) - 0.5)*(Lx)
values = func(points[:,0], points[:,1])
values = np.array([values])
positions = np.concatenate((points, values.T), axis = 1)
z = positions[:,2]
r = positions[:,0:2]
[udirect_slow] = compute_u_direct_slow(Lx, Ly, Ndirect, r, z)'''
    
    # timeit.repeat statement
    times = timeit.repeat(setup = SETUP_CODE,
                          stmt = TEST_CODE,
                          repeat = 3,
                          number = 100)
    print('Slow time: {}'.format(min(times)))

# Compute the fast time
def fast_time():
    SETUP_CODE = '''
from __main__ import compute_u_direct_fast, func
import numpy as np
'''

    TEST_CODE = '''
Nx = 200
Ny = 200
Lx = 200
Ly = 200
npoints = 2000
Ndirect = 21

rng = np.random.default_rng()
points = (rng.random((npoints, 2)) - 0.5)*(Lx)
values = func(points[:,0], points[:,1])
values = np.array([values])
positions = np.concatenate((points, values.T), axis = 1)
z = positions[:,2]
r = positions[:,0:2]
[udirect_fast] = compute_u_direct_fast(Lx, Ly, Ndirect, r, z)'''
    
    # timeit.repeat statement
    times = timeit.repeat(setup = SETUP_CODE,
                          stmt = TEST_CODE,
                          repeat = 3,
                          number = 100)
    print('Fast time: {}'.format(min(times)))
    


if __name__ == "__main__":
    #fft_time()
    #slow_time()
    #fast_time()

    # Set up the system
    Nx = 200
    Ny = 200
    Lx = 200
    Ly = 200
    Nxgrid = Nx * 1j
    Nygrid = Ny * 1j
    npoints = 2000
    Ndirect = 11

    # Set the direct calc cutoff
    qcutoff_direct = 2.0*np.pi/Lx
    
    # Create some random points
    rng = np.random.default_rng()
    points = (rng.random((npoints, 2)) - 0.5)*(Lx)
    values = func(points[:,0], points[:,1])
    values = np.array([values])
    
    # Concatenate together to make something like a position statement
    positions = np.concatenate((points, values.T), axis = 1)
    
    z = positions[:,2]
    r = positions[:,0:2]
    
    ################
    # Interpolation method on a grid
    ################
    
    # Save the information
    [u_fft, qcutoff_fft] = compute_u_fft(Lx, Ly, Nx, Ny, r, z)
    
    ################
    # Original direct (slow) version
    ################
    [udirect_slow] = compute_u_direct_slow(Lx, Ly, Ndirect, r, z)
    
    ################
    # Fast direct version?
    ################
    [udirect_fast] = compute_u_direct_fast(Lx, Ly, Ndirect, r, z)
    
    ################
    # Plot everythig to compare
    ################
    fig, axarr = plt.subplots(1, 3)
    axarr[0].imshow(np.abs(u_fft.T), origin = 'lower')
    axarr[1].imshow(np.abs(udirect_slow.T), origin = 'lower')
    axarr[2].imshow(np.abs(udirect_fast.T), origin = 'lower')
    
    plt.show()

    # Up to here, everything looks good, except for the normalization factor. Need to figure
    # out the difference in the scaling of the real and imaginary parts

    # Find the maximum value in abs space, and then use that to draw a 1D line across each one to compare
    idx_fft     = np.unravel_index(np.argmax(u_fft, axis=None), u_fft.shape)
    idx_slow    = np.unravel_index(np.argmax(udirect_slow, axis=None), udirect_slow.shape)
    idx_fast    = np.unravel_index(np.argmax(udirect_fast, axis=None), udirect_fast.shape)


    print(f"Bailing early to figure out scaling argument")
    sys.exit(1)

    # Save everything to compare radial averaging
    uq_2d_fft = u_fft
    uq_2d_slow = udirect_slow
    uq_2d_fast = udirect_fast

    # Chose a deltaq
    deltaq = 0.0375

    [radii_fft, intensity_fft] = radial_average(uq_2d_fft, deltaq, qcutoff_fft)
    [radii_slow, intensity_slow] = radial_average(uq_2d_slow, deltaq, qcutoff_direct)
    [radii_fast, intensity_fast] = radial_average(uq_2d_fast, deltaq, qcutoff_direct)

    su_fft  = np.square(intensity_fft)*len(z)
    su_slow = np.square(intensity_slow)*len(z)
    su_fast = np.square(intensity_fast)*len(z)

    fig, axarr = plt.subplots(1, 1)

    axarr.scatter(radii_fft, su_fft, color = 'b', marker = '+', linewidth = 1)
    axarr.scatter(radii_slow, su_slow, color = 'r', marker = 's', s = 80, facecolors = 'none')
    axarr.scatter(radii_fast, su_fast, color = 'm', marker = 'o', s = 80, facecolors = 'none')

    axarr.set_yscale('log')
    axarr.set_xscale('log')

    fig.tight_layout()
    plt.show()
