#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import scipy
import timeit

# Defined function for use
def func(x, y):
    return np.cos(4*np.pi*x/100) * np.sin(8*np.pi*y/100)

# Compute the FFT version
def compute_u_fft(Lx, Ly, Nx, Ny, r, z):
    xpixel = Lx/Nx
    ypixel = Ly/Ny
    Nxgrid = Nx * 1j
    Nygrid = Ny * 1j
    grid_x, grid_y = np.mgrid[-Lx/2:Lx/2:Nxgrid, -Ly/2:Ly/2:Nygrid]
    from scipy.interpolate import griddata
    grid_z = griddata(r, z, (grid_x, grid_y), method = 'cubic', fill_value = np.mean(z))

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
    fft_time()
    #slow_time()
    fast_time()


## Set up the system
#Nx = 200
#Ny = 200
#Lx = 200
#Ly = 200
#Nxgrid = Nx * 1j
#Nygrid = Ny * 1j
#npoints = 2000
#Ndirect = 1
#
## Create some random points
#rng = np.random.default_rng()
#points = (rng.random((npoints, 2)) - 0.5)*(Lx)
#values = func(points[:,0], points[:,1])
#values = np.array([values])
#
## Concatenate together to make something like a position statement
#positions = np.concatenate((points, values.T), axis = 1)
#
#z = positions[:,2]
#r = positions[:,0:2]
#
#################
## Interpolation method on a grid
#################
#
## Save the information
#[u_fft, qcutoff_fft] = compute_u_fft(Lx, Ly, Nx, Ny, r, z)
#
#################
## Original direct (slow) version
#################
#[udirect_slow] = compute_u_direct_slow(Ndirect, r, z)
#
#################
## Fast direct version?
#################
#[udirect_fast] = compute_u_direct_fast(Ndirect, r, z)
#
#################
## Plot everythig to compare
#################
#plt.subplot(131)
#plt.imshow(np.abs(u_fft.T), origin = 'lower')
#plt.title('Cubic FFT')
#plt.subplot(132)
#plt.imshow(np.abs(udirect_slow.T), origin = 'lower')
#plt.title('Direct (slow)')
#plt.subplot(133)
#plt.imshow(np.abs(udirect_fast.T), origin = 'lower')
#plt.title('Direct (fast)')
#
#plt.show()
