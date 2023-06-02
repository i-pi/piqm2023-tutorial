import re
import numpy as np


# Functions for evaluating the mean energy for the harmonic potential, analytically
def getZk(k, bhw, dim):
    return np.power((np.exp(0.5 * k * bhw) / (np.exp(k * bhw) - 1)), dim)


def getdZk(k, bhw, dim):
    return -0.5 * k * dim * getZk(k, bhw, dim) * (1 + np.exp(-k * bhw)) / (1 - np.exp(-k * bhw))


# Mean energy of n non-interacting bosons (Alg. 4.7 from Krauth, p. 231)
def get_harmonic_energy(n, bhw, dim, is_bosonic=True):
    if not is_bosonic:
        return -n * getdZk(1, bhw, dim) / getZk(1, bhw, dim)

    z_arr = np.zeros(n + 1)
    dz_arr = np.zeros(n + 1)

    z_arr[0] = 1.0

    for m in range(1, n + 1):
        sig_z = 0.0
        sig_dz = 0.0

        for j in range(m, 0, -1):
            sig_z += getZk(j, bhw, dim) * z_arr[m - j]
            sig_dz += getdZk(j, bhw, dim) * z_arr[m - j] + getZk(j, bhw, dim) * dz_arr[m - j]

        z_arr[m] = sig_z / m
        dz_arr[m] = sig_dz / m

    return -dz_arr[n] / z_arr[n]


def read_ipi_output(filename):
    """ Reads an i-PI output file and returns a dictionary with the properties in a tidy order. """

    f = open(filename, "r")

    regex = re.compile(".*column *([0-9]*) *--> ([^ {]*)")

    fields = []
    cols = []
    for line in f:
        if line[0] == "#":
            match = regex.match(line)
            if match is None:
                print("Malformed comment line: ", line)
                continue  # TODO: was error
            fields.append(match.group(2))
            cols.append(slice(int(match.group(1)) - 1, int(match.group(1))))
        else:
            break  # done with header
    f.close()

    columns = {}
    raw = np.loadtxt(filename)
    for i, c in enumerate(fields):
        columns[c] = raw[:, cols[i]].T
        if columns[c].shape[0] == 1:
            columns[c].shape = columns[c].shape[1]
    return columns


def get_mean_energies(filename, skip_steps=100):
    o = read_ipi_output(filename + "data.out")
    avg_kinetic = -o['virial_fq'][skip_steps:]
    avg_kinetic_cv = o['kinetic_cv'][skip_steps:]
    avg_potential = o['potential'][skip_steps:]
    return {
        'virial': avg_kinetic,
        'centroid_virial': avg_kinetic_cv,
        'potential': avg_potential
    }


def analytical_energy(temp=17.4, sys_type='dist'):
    dim = 3  # Dimension of the physical system
    spring_constant = 1.21647924E-8  # Spring constant
    mass = 1.0  # Mass of the particle

    # hbar and the Boltzmann constant are assumed to be 1
    hbar = 1
    kB = 1
    omega = np.sqrt(spring_constant / mass)

    # See units.py of i-pi for conversion to atomic units
    beta = (1 / (kB * temp * 3.1668152e-06))

    bhw = hbar * beta * omega

    if sys_type == 'bosonic':
        return get_harmonic_energy(3, bhw, dim, True) * hbar * omega
    elif sys_type == 'mixed':
        return (get_harmonic_energy(2, bhw, dim, True) + get_harmonic_energy(1, bhw, dim, False)) * hbar * omega
    elif sys_type == 'fermionic':
        # Analytical result for three fermions in a 2D harmonic trap (NVT)
        exp = np.exp(bhw)
        exp2 = np.exp(2 * bhw)
        exp3 = np.exp(3 * bhw)
        exp4 = np.exp(4 * bhw)
        exp5 = np.exp(5 * bhw)
        exp6 = np.exp(6 * bhw)
        num =  5 * exp6 + 31 * exp5 + 47 * exp4 + 50 * exp3 + 47 * exp2 + 31 * exp + 5
        denom = (exp - 1) * (exp + 1) * (exp2 + exp + 1) * (exp2 + 4 * exp + 1)

        return hbar * omega * num / denom

    # Unless specified otherwise, assume distinguishable particles.
    return get_harmonic_energy(3, bhw, dim, False) * hbar * omega