#!/bin/bash

# Read mass values after time-step 10
m_0=$(awk 'NR == 11 {print $1}' Moments_Test0.dc)
m=$(awk 'NR == 11 {print $1}' Data/Moments_nu0.05A0.2k0.5Nx16Lx12.5664Nv16Lv5.25SpectralN8dt0.01nT10_Test0.dc)

if [ $(echo "$m - $m_0 > 0.000002" | bc) == 1 ] || [ $(echo "$m_0 - $m > 0.000002" | bc) == 1 ]; then
	echo "Error in mass value!"
	exit 1
fi

# Read first component of velocity values after time-step 10
u1_0_in=$(awk 'NR == 11 {print $2}' Moments_Test0.dc)
u1_0=$(echo "$u1_0_in" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
u1_in=$(awk 'NR == 11 {print $2}' Data/Moments_nu0.05A0.2k0.5Nx16Lx12.5664Nv16Lv5.25SpectralN8dt0.01nT10_Test0.dc)
u1=$(echo "$u1_in" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')

if [ $(echo "$u1 - $u1_0 > 10^-10" | bc -l) == 1 ] || [ $(echo "$u1_0 - $u1 > 10^-10" | bc -l) == 1 ]; then
	echo "Error in first component of velocity!"
	exit 1
fi

# Read second component of velocity values after time-step 10
u2_0_in=$(awk 'NR == 11 {print $3}' Moments_Test0.dc)
u2_0=$(echo "$u2_0_in" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
u2_in=$(awk 'NR == 11 {print $3}' Data/Moments_nu0.05A0.2k0.5Nx16Lx12.5664Nv16Lv5.25SpectralN8dt0.01nT10_Test0.dc)
u2=$(echo "$u2_in" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')

if [ $(echo "$u2 - $u2_0 > 10^-10" | bc -l) == 1 ] || [ $(echo "$u2_0 - $u2 > 10^-10" | bc -l) == 1 ]; then
	echo "Error in second component of velocity!"
	exit 1
fi

# Read third component of velocity values after time-step 10
u3_0_in=$(awk 'NR == 11 {print $4}' Moments_Test0.dc)
u3_0=$(echo "$u3_0_in" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
u3_in=$(awk 'NR == 11 {print $4}' Data/Moments_nu0.05A0.2k0.5Nx16Lx12.5664Nv16Lv5.25SpectralN8dt0.01nT10_Test0.dc)
u3=$(echo "$u3_in" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')

if [ $(echo "$u3 - $u3_0 > 10^-10" | bc -l) == 1 ] || [ $(echo "$u3_0 - $u3 > 10^-10" | bc -l) == 1 ]; then
	echo "Error in third component of velocity!"
	exit 1
fi

# Read total energy values after time-step 10
T_0=$(awk 'NR == 11 {print $9}' Moments_Test0.dc)
T=$(awk 'NR == 11 {print $9}' Data/Moments_nu0.05A0.2k0.5Nx16Lx12.5664Nv16Lv5.25SpectralN8dt0.01nT10_Test0.dc)
if [ $(echo "$T - $T_0 > 0.00003" | bc) == 1 ] || [ $(echo "$T_0 - $T > 0.00003" | bc) == 1 ]; then
	echo "Error in total energy!"
	exit 1
fi

exit 0
