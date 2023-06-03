import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import MDAnalysis as mda
from pathlib import Path
from analysis import *


def position_tensor(beads, N, filesPath):
    """
    READ FILES IN XYZ from i-PI
    notice a and b are just dummy arrays.
    Returns
    -------
    Position tensor: R[steps, beads, N, x, y, z].
    """
    conv_angst_to_bohr = 0.529177  # 1 bohr = 0.529177 angstrom

    path = Path(filesPath)

    # Read XYZ files and fetch all the coordinates at all the timesteps
    xyz_filenames = list(path.rglob("*.xyz"))

    xyz_class = mda.coordinates.XYZ.XYZReader(xyz_filenames[0])

    steps = xyz_class.n_frames

    # coords[timestep][bead_idx][ptcl_idx][axis]
    coords = np.zeros((steps, beads, N, 3))

    for bead_idx, bead_file in enumerate(xyz_filenames):
        xyz_file = mda.coordinates.XYZ.XYZReader(bead_file)

        for idx, timestep in enumerate(xyz_file):
            coords[timestep.frame, bead_idx, :, :] = timestep.positions

    return coords * conv_angst_to_bohr, steps


class ExchangePotential():
    def __init__(
        self, N, q, nbeads, filesPath, spring_freq_squared, bead_mass
    ):
        self._N = N
        self._P = nbeads
        self._q = q
        self._spring_freq_squared = spring_freq_squared
        self._particle_mass = bead_mass
        self._filesPath = filesPath

        # self._bead_diff_intra[j] = [r^{j+1}_0 - r^{j}_0, ..., r^{j+1}_{N-1} - r^{j}_{N-1}]
        self._bead_diff_intra = np.diff(self._q, axis=0)
        # self._bead_dist_inter_first_last_bead[l][m] = r^0_{l} - r^{P-1}_{m}
        self._bead_diff_inter_first_last_bead = (
            self._q[0, :, np.newaxis, :] - self._q[self._P - 1, np.newaxis, :, :]
        )

        # cycle energies:
        # self._E_from_to[u, u] is the ring polymer energy of the cycle on particle indices u,...,v
        self._E_from_to = self._evaluate_cycle_energies()

    def _spring_potential_prefix(self):
        """
        Helper function: the term for the spring constant as used in spring potential expressions
        """
        return 0.5 * self._particle_mass * self._spring_freq_squared


    def _evaluate_cycle_energies(self):
        """
        Evaluate all the cycle energies, as outlined in Eqs. 5-7 of the paper.
        Returns an upper-triangular matrix, Emks[u,v] is the ring polymer energy of the cycle on u,...,v.
        """
        # using column-major (Fortran order) because uses of the array are mostly within the same column
        Emks = np.zeros((self._N, self._N), dtype=float, order="F")

        intra_spring_energies = np.sum(self._bead_diff_intra ** 2, axis=(0, -1))
        spring_energy_first_last_bead_array = np.sum(
            self._bead_diff_inter_first_last_bead ** 2, axis=-1
        )

        # for m in range(self._N):
        #     Emks[m][m] = self._spring_potential_prefix() * \
        #                     (intra_spring_energies[m] + spring_energy_first_last_bead_array[m, m])
        Emks[np.diag_indices_from(Emks)] = self._spring_potential_prefix() * (
            intra_spring_energies + np.diagonal(spring_energy_first_last_bead_array)
        )

        for s in range(self._N - 1 - 1, -1, -1):
            # for m in range(s + 1, self._N):
            #     Emks[s][m] = Emks[s + 1][m] + self._spring_potential_prefix() * (
            #             - spring_energy_first_last_bead_array[s + 1, m]
            #             + intra_spring_energies[s]
            #             + spring_energy_first_last_bead_array[s + 1, s]
            #             + spring_energy_first_last_bead_array[s, m])
            Emks[s, (s + 1) :] = Emks[
                s + 1, (s + 1) :
            ] + self._spring_potential_prefix() * (
                -spring_energy_first_last_bead_array[s + 1, (s + 1) :]
                + intra_spring_energies[s]
                + spring_energy_first_last_bead_array[s + 1, s]
                + spring_energy_first_last_bead_array[s, (s + 1) :]
            )

        return Emks


def total_fermionic_energy(path, sign_arr, skip_steps, skip_freq=1):
    o = read_ipi_output(path + "/data.out")
    ke_list = -o['virial_fq'][skip_steps:]
    pe_list = o['potential'][skip_steps:]

    # Adjust skip_freq depending on thermo and dump of system
    tot_est = pe_list[::skip_freq] + ke_list[::skip_freq]

    e_s_num_h = tot_est * sign_arr
    e_s_num_h = np.mean(e_s_num_h)
    s_denom = np.mean(sign_arr)

    avg_etot_f_h = e_s_num_h / s_denom

    return avg_etot_f_h


def get_enk(Enk_tensor, t, N, k):
    if N == 0:
        return 0

    return Enk_tensor[t][N - k][N - 1]


def get_W(Enk_tensor, t, beta, N, fermionic=True):
    xi = -1 if fermionic else 1
    w_arr = np.zeros(N + 1)
    w_arr[0] = 1  # Initial condition for the recurrence relation

    for m in range(1, N + 1):
        w_val_tot = 0
        for k in range(m, 0, -1):
            w_val = xi**(k - 1)
            w_val *= np.exp(-beta * get_enk(Enk_tensor, t, m, k))
            w_val *= w_arr[m - k]
            w_val_tot += w_val

        w_arr[m] = w_val_tot / m

    return w_arr[N]


def get_reweighted_energy(path, N=3, P=12, T=30, skip_steps=500):
    hbar = 1
    # hw = 1.1e-4  # In Hartree (=0.003 eV)
    m = 1  # Atomic mass units
    kB = 3.16683 * 10 ** (-6)  # Boltzmann in Hartree/K
    # kB = 8.61733e-5  # ev/K
    beta = 1 / (kB * T)
    w2 = P / (beta ** 2 * hbar ** 2)

    R, steps = position_tensor(P, N, path)
    enk_array = np.array([])

    for t in range(int(steps)):
        enk = ExchangePotential(N, R[t, :, :, :], P, path, w2, m)
        enk_array = np.append(enk_array, enk._E_from_to.flatten())

    enk_array = enk_array.reshape(steps, N, N)

    sign_array = []
    for t in range(steps):
        sign_array.append(get_W(enk_array, t, beta, N, True) / get_W(enk_array, t, beta, N, False))

    sign_array = sign_array[skip_steps:]

    avg_etot_f_h = total_fermionic_energy(path, np.array(sign_array), skip_steps, 1)

    return np.mean(sign_array), avg_etot_f_h