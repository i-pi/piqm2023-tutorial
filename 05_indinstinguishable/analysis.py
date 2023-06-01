import re
import sys

import numpy as np
import statistics
import argparse
import pathlib
import itertools

#import harmonic_analytical

T = 17.8  # Kelvin
nparticles = 3

dim = 3  # Dimension of the physical system
spring_constant = 1.21647924E-8  # Spring constant
mass = 1.0  # Mass of the particle

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


def analyze_mean_energy(infile):
    o = read_ipi_output(infile)
    num_points_to_drop = 100
    avg_kinetic = -statistics.mean(o['virial_fq'][num_points_to_drop:])
    avg_potential = statistics.mean(o['potential'][num_points_to_drop:])
    return avg_kinetic + avg_potential


def harmonic_analytical_energy(nparticles, temp, bosonic=True):
    # hbar and the Boltzmann constant are assumed to be 1

    omega = np.sqrt(spring_constant / mass)

    # See units.py of i-pi for conversion to atomic units
    beta = (1 / (temp * 3.1668152e-06))

    bhw = beta * omega

    return get_harmonic_energy(nparticles, bhw, dim, bosonic)


def get_average_energy(infiles, particles, temp):
    avg_energy = analyze_mean_energy(infiles[0])
    analytical_energy = harmonic_analytical_energy(particles, temp)

    print("Avg energy:", avg_energy)
    print("Analytical:", analytical_energy)
    return avg_energy


def main():
    parser = argparse.ArgumentParser(description='Analyze i-Pi simulation results')

    parser.add_argument('infile', type=str, nargs='*', help='Path to i-Pi data.out output file')
    parser.add_argument('-t', type=float, nargs=1, help='Temperature (in kelvin)')

    # Think about how to deal with case of mixed particles and/or fermions
    parser.add_argument('-n', type=int, nargs=1, help='Number of particles')

    args = parser.parse_args()

    # analyze_mean_energy(args.infile[0])
    get_average_energy(args.infile, args.n[0], args.t[0])


if __name__ == "__main__":
    main()