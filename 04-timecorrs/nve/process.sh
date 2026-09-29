#!/bin/bash
#
# Computes the velocity autocorrelation function (VACF) and its Fourier
# transform (vibrational power spectrum) from the NVE trajectory.
#
# Reads simulation.vel_0.xyz and uses i-pi-getacf to:
#   - compute the VACF out to a maximum lag of 1024 steps (-mlag)
#   - zero-pad the VACF by 3072 points before the FFT to increase the
#     frequency resolution of the spectrum (-ftpad)
#   - apply a cosine-Hanning window to the VACF before the FFT, to
#     reduce spectral leakage from the finite trajectory length (-ftwin)
#   - assume a trajectory frame spacing of 1.0 fs (-dt)
#
# Outputs (prefix "nve"):
#   nve_acf.data  : columns [time, VACF, VACF_error]
#   nve_facf.data : columns [angular frequency, power spectrum, power spectrum error]

i-pi-getacf -ifile simulation.vel_0.xyz -mlag 1024 -ftpad 3072 -ftwin cosine-hanning -dt "1.0 femtosecond" -oprefix nve

gnuplot -persist <<'EOF'
set terminal qt size 600,600 persist
set multiplot layout 2,1

# unit conversions from atomic units 
au2fs = 1.0 / 41.341373      # atomic time unit -> femtosecond
au2wn = 1.0 / 4.5563353e-06  # atomic angular frequency unit -> wavenumber (cm^-1)

set title "Velocity autocorrelation function"
set xlabel "Time (fs)"
set ylabel "VACF (a.u.)"
set xrange [0:1000]
plot "nve_acf.data" using ($1*au2fs):2 with lines notitle

set title "Vibrational power spectrum"
set xlabel "Wavenumber (cm^{-1})"
set ylabel "Power spectrum (a.u.)"
set xrange [0:4500]
plot "nve_facf.data" using ($1*au2wn):2 with lines notitle

unset multiplot
EOF
