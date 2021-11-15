#!/usr/bin/env python3

import fresnel
import freud
import gsd.hoomd
import sys

import numpy as np
import matplotlib.pyplot as plt

# Create a fake signal
def func(x, y):
    return np.cos(2*np.pi*x/100) * np.sin(4*np.pi*y/100)

xsize = 100
ysize = 100
grid_x, grid_y = np.mgrid[-xsize/2:xsize/2:200j,-ysize/2:ysize/2:200j]

# Create some random points
rng = np.random.default_rng()
points = (rng.random((1000, 2)) - 0.5)*(xsize)
values = func(points[:,0], points[:,1])

# Create the grid values
from scipy.interpolate import griddata
grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')

plt.subplot(221)
plt.imshow(func(grid_x, grid_y).T, extent=(-xsize/2, xsize/2, -ysize/2, ysize/2), origin='lower')
plt.plot(points[:,0], points[:,1], 'k.', ms=1)
plt.title('Original')
plt.subplot(222)
plt.imshow(grid_z0.T, extent=(0,1,0,1), origin='lower')
plt.title('Nearest')
plt.subplot(223)
plt.imshow(grid_z1.T, extent=(0,1,0,1), origin='lower')
plt.title('Linear')
plt.subplot(224)
plt.imshow(grid_z2.T, extent=(0,1,0,1), origin='lower')
plt.title('Cubic')
plt.gcf().set_size_inches(6, 6)
plt.show()

# Now take the 2D fourier transform of this
# FFT based on the nearest neighbor
u = np.fft.fft2(grid_z0, norm='forward')
ushift = np.fft.fftshift(u)
sfreq = np.fft.fftfreq(xsize)
sfreqshift = np.fft.fftshift(sfreq)
cen_x = np.int32(ushift.shape[1]/2)
cen_y = np.int32(ushift.shape[0]/2)

# FFT of the cubic interpolation
grid_z2[np.isnan(grid_z2)] = np.mean(values)
ucubic = np.fft.fft2(grid_z2, norm='forward')
ucubicshift = np.fft.fftshift(ucubic)
ccen_x = np.int32(ucubic.shape[1]/2)
ccen_y = np.int32(ucubic.shape[0]/2)

# Directly do the transform and see how it compares
Nx = 20
deltaq = 1.0
udirect = np.zeros((2*Nx, 2*Nx), dtype=np.complex128)
for n in np.arange(-Nx,Nx,1):
    for m in np.arange(-Nx,Nx,1):
        idx = n + Nx
        jdx = m + Nx
        
        q = 2.0*np.pi*np.array([n/xsize, m/ysize])
        qnorm = np.linalg.norm(q)
        #print(f"({n}, {m})({idx}, {jdx}) -> ({q[0]}, {q[1]}) -> qnorm: {qnorm}")
        for k,z1k in enumerate(values):
            r1k = np.array([points[k,0], points[k,1]])
            val = z1k * np.exp(-1j*np.dot(q, r1k))
            udirect[idx,jdx] += val

# Normalize the direct version
udirect = udirect / len(values)
dcen_x = np.int32(udirect.shape[1]/2)
dcen_y = np.int32(udirect.shape[0]/2)

plt.subplot(131)
plt.imshow(np.abs(ushift.T), origin='lower')
plt.title('Grid FFT')
plt.subplot(132)
plt.imshow(np.abs(ucubicshift.T), origin='lower')
plt.title('Cubic FFT')
plt.subplot(133)
plt.imshow(np.abs(udirect.T), origin='lower')
plt.title('Direct FT')
plt.show()

# Try to do some radial averaging
# Fourier grid spacing one
X,Y = np.meshgrid(np.arange(ushift.shape[1]) - cen_x, np.arange(ushift.shape[0]) - cen_y)
R = np.sqrt(X**2 + Y**2)
rad = np.arange(1, np.max(R), 1)
intensity = np.zeros(len(rad))
intensity_cubic = np.zeros(len(rad))
index = 0
bin_size = 1.0
for i in rad:
    mask = (np.greater(R, i - bin_size) & np.less(R, i + bin_size))
    values = np.abs(ushift[mask])
    values_cubic = np.abs(ucubicshift[mask])
    intensity[index] = np.mean(values)
    intensity_cubic[index] = np.mean(values_cubic)
    index += 1

# Direct calculation
X,Y = np.meshgrid(np.arange(udirect.shape[1]) - dcen_x, np.arange(udirect.shape[0]) - dcen_y)
R = np.sqrt(X**2 + Y**2)
rad_direct = np.arange(1, np.max(R), 1)
intensity_direct = np.zeros(len(rad_direct))
index = 0
bin_size = 1.0
for i in rad_direct:
    mask = (np.greater(R, i - bin_size) & np.less(R, i + bin_size))
    values = np.abs(udirect[mask])
    intensity_direct[index] = np.mean(values)
    index += 1

plt.subplot(111)
plt.loglog(rad, np.square(intensity), 'r+')
plt.loglog(rad, np.square(intensity_cubic), 'r^')
plt.loglog(rad_direct, np.square(intensity_direct), 'bo')

plt.show()



