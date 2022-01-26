# Common helper functions for membranes and AH domains

import random
import numpy as np

# Generate a random unit vector
def generate_random_unit_vector(ndim):
    w = np.float64(1.0)
    v = np.zeros(3, dtype = np.float64)
    if ndim == 3:
        z = np.float64(2.0 * random.random() - 1.0)
        w = np.sqrt(1.0 - z*z)
        v[2] = z
    t = 2.0*np.pi*random.random()
    v[0] = w*np.cos(t)
    v[1] = w*np.sin(t)
    return v

# Check for an sphere overlap
def sphere_overlap(sphere_diameter, r1, r2):
    distmag = np.linalg.norm(r1 - r2)
    return distmag < sphere_diameter

# Define a moving average
def moving_average(a, n=2):
    r""" Moving average
    """
    ret = np.cumsum(a, dtype=np.float32)
    ret[n:] = ret[n:] - ret[:-n]
    return np.divide(ret[n-1:],n)

# Define a radial average
def radial_average(uq, bin_size, pixel_size):
    cen_x = np.int32(uq.shape[1]/2)
    cen_y = np.int32(uq.shape[0]/2)
    X, Y = np.meshgrid(np.arange(uq.shape[1]) - cen_x, np.arange(uq.shape[0]) - cen_y)
    R = np.sqrt(X**2 + Y**2)
    # Rreal has the information about the actual location of points
    Rreal = pixel_size * R

    # Set Rmax so that we don't go outside the circle averaging
    Rmax = Rreal[0, cen_x]

    # Calculate the number of bins
    nbins = np.int32(np.floor(Rmax / bin_size))
    bin_edges = np.arange(0, nbins+1)*bin_size
    bin_mids = (bin_edges[0:-1] + bin_edges[1:])/2

    intensity = np.zeros(len(bin_mids))
    # Loop over the bins and count up what we need
    for ibin in np.arange(nbins):
        mask = (np.greater_equal(Rreal, bin_edges[ibin])) & (np.less(Rreal, bin_edges[ibin+1]))
        values = np.abs(uq[mask])
        intensity[ibin] = np.mean(values)

    return [bin_mids, intensity]
