#!/usr/bin/env python3

import fresnel
import freud
import gsd.hoomd
import sys

import numpy as np
import matplotlib.pyplot as plt

import scipy

# Create a fake signal
def func(x, y):
    return np.cos(4*np.pi*x/100) * np.sin(8*np.pi*y/100)

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
        intensity[ibin] = np.mean(values)

    return [bin_mids, intensity]

# Create some fake points
Lx = 200
Ly = 200
Nx = 200
Ny = 200
Nxgrid = 200j
Nygrid = 200j

grid_x, grid_y = np.mgrid[-Lx/2:Lx/2:Nxgrid,-Ly/2:Ly/2:Nygrid]
xpixel = Lx/Nx
ypixel = Ly/Nx

# Create some random points
rng = np.random.default_rng()
points = (rng.random((2000, 2)) - 0.5)*(Lx)
values = func(points[:,0], points[:,1])

# Create the grid values
from scipy.interpolate import griddata
grid_z1 = griddata(points, values, (grid_x, grid_y), method='cubic', fill_value = np.mean(values))

plt.subplot(121)
plt.imshow(func(grid_x, grid_y).T, extent=(-Lx/2, Lx/2, -Ly/2, Ly/2), origin='lower')
plt.plot(points[:,0], points[:,1], 'k.', ms=1)
plt.title('Original')
plt.subplot(122)
plt.imshow(grid_z1.T, extent=(-Lx/2, Lx/2, -Ly/2, Ly/2), origin='lower')
plt.title('Cubic')
plt.gcf().set_size_inches(6, 6)
plt.show()

# Now take the 2D fourier transform of this
# FFT based on the nearest neighbor
ucubic = np.fft.fft2(grid_z1, norm='forward')
ucubicshift = np.fft.fftshift(ucubic)
freqx = np.fft.fftshift(np.fft.fftfreq(ucubicshift.shape[1], xpixel))
freqy = np.fft.fftshift(np.fft.fftfreq(ucubicshift.shape[0], ypixel))
ccen_x = np.int32(ucubic.shape[1]/2)
ccen_y = np.int32(ucubic.shape[0]/2)

# Directly do the transform and see how it compares
Ndirect = 21
udirect = np.zeros((2*Ndirect+1, 2*Ndirect+1), dtype=np.complex128)
qdirectx = np.zeros(2*Ndirect+1, dtype=np.float64)
for n in np.arange(-Ndirect, Ndirect+1, 1):
    for m in np.arange(-Ndirect, Ndirect+1, 1):
        idx = n + Ndirect
        jdx = m + Ndirect
        
        q = 2.0*np.pi*np.array([n/Lx, m/Lx])
        qnorm = np.linalg.norm(q)
        qdirectx[idx] = q[0]
        print(f"({n}, {m})({idx}, {jdx}) -> ({q[0]}, {q[1]}) -> qnorm: {qnorm}")
        for k,z1k in enumerate(values):
            r1k = np.array([points[k,0], points[k,1]])
            val = z1k * np.exp(-1j*np.dot(q, r1k))
            udirect[idx,jdx] += val

# Normalize the direct version
udirect = udirect / len(values)

plt.subplot(121)
plt.imshow(np.abs(ucubicshift.T), origin='lower')
plt.title('Cubic FFT')
plt.subplot(133)
plt.imshow(np.abs(udirect.T), origin='lower')
plt.title('Direct FT')
plt.show()

# What do we think the frequencies are?
print(f"FFT frequencies")
print(f"{freqx}")
print(f"Direct wave number q")
print(f"{qdirectx}")

freq_direct_pixel = qdirectx[1] - qdirectx[0]
freq_fft_pixel = (freqx[1] - freqx[0])*(2.0*np.pi)

print(f"Frequency direct pixel: {freq_direct_pixel}")
print(f"Frequency FFT pixel:    {freq_fft_pixel}")

# Now do the radial averaging
deltaq = 0.1
[radii_direct, intensity_direct] = radial_average(udirect, deltaq, freq_direct_pixel)
[radii_fft, intensity_fft] = radial_average(ucubicshift, deltaq, freq_fft_pixel)

print(f"Direct:")
print(f"  radii_direct:     {radii_direct}")
print(f"  intensity_direct: {intensity_direct}")
print(f"FFT:")
print(f"  radii_fft:        {radii_fft}")
print(f"  intensity_fft:    {intensity_fft}")

# Try to plot this mess
plt.subplot(111)
plt.loglog(radii_direct, np.square(intensity_direct), 'bo')
plt.loglog(radii_fft, np.square(intensity_fft), 'r+')
plt.axvline(x =  2*np.pi/Lx,  ymin = 0, ymax = 1, color = 'red', linestyle = '-.')
plt.axvline(x =  4*np.pi/100, ymin = 0, ymax = 1, color = 'black', linestyle = '--')
plt.axvline(x =  8*np.pi/100, ymin = 0, ymax = 1, color = 'black', linestyle = ':')
plt.show()


## Try to do some radial averaging
## Fourier grid spacing one
#X,Y = np.meshgrid(np.arange(ushift.shape[1]) - cen_x, np.arange(ushift.shape[0]) - cen_y)
#R = np.sqrt(X**2 + Y**2)
#rad = np.arange(1, np.max(R), 1)
#intensity = np.zeros(len(rad))
#intensity_cubic = np.zeros(len(rad))
#index = 0
#bin_size = 1.0
#for i in rad:
#    mask = (np.greater(R, i - bin_size) & np.less(R, i + bin_size))
#    values = np.abs(ushift[mask])
#    values_cubic = np.abs(ucubicshift[mask])
#    intensity[index] = np.mean(values)
#    intensity_cubic[index] = np.mean(values_cubic)
#    index += 1
#
## Direct calculation
#X,Y = np.meshgrid(np.arange(udirect.shape[1]) - dcen_x, np.arange(udirect.shape[0]) - dcen_y)
#R = np.sqrt(X**2 + Y**2)
#rad_direct = np.arange(1, np.max(R), 1)
#intensity_direct = np.zeros(len(rad_direct))
#index = 0
#bin_size = 1.0
#for i in rad_direct:
#    mask = (np.greater(R, i - bin_size) & np.less(R, i + bin_size))
#    values = np.abs(udirect[mask])
#    intensity_direct[index] = np.mean(values)
#    index += 1
#
#plt.subplot(111)
#plt.loglog(rad, np.square(intensity), 'r+')
#plt.loglog(rad, np.square(intensity_cubic), 'r^')
#plt.loglog(rad_direct, np.square(intensity_direct), 'bo')
#
#plt.show()
#
#
#
