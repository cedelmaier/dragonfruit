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

