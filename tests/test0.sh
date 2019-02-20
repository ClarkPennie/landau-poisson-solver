#!/bin/bash

#==============================================
# Test for the Landau Damping ICs
#==============================================
echo "Testing Landau Damping ICs..."

# Set the name of the files with the moment values to be compared
moment_filename_expected=Moments_Test0.dc
moment_filename_test=Data/Moments_nu0.05A0.2k0.5Nx16Lx12.5664Nv16Lv5.25SpectralN8dt0.01nT5_Test0.dc

# Remove the file containing the moments to be produced in case it already exists
rm $moment_filename_test

echo "Running code for 5 timesteps..." 
mv LPsolver-input-test0.txt LPsolver-input.txt
../source/solver 
mv LPsolver-input.txt LPsolver-input-test0.txt

echo "Checking values of moments are as expected..."
sh moment_differ.sh $moment_filename_expected $moment_filename_test
exit_val=$?

if [ $exit_val -eq 1 ]; then
	exit $exit_val
fi

#==============================================
# Test for the non-uniform doping ICs
#==============================================
echo "Testing non-uniform doping ICs..."

# Set the name of the files with the moment values to be compared
moment_filename_expected=Moments_Test1.dc
moment_filename_test=Data/Moments_nu0.05A0k0.5Nx16Lx12.5664Nv16Lv5.25SpectralN8dt0.01nT5_Test1.dc

# Remove the file containing the moments to be produced in case it already exists
rm $moment_filename_test

# This version does not conserve all moments, so set new thresholds
m_thresh=0.000002
u1_thresh=0.01
u2_thresh=2*10^-7
u3_thresh=2*10^-7
T_thresh=0.01

echo "Running code for 5 timesteps..." 
mv LPsolver-input-test1.txt LPsolver-input.txt
../source/solver 
mv LPsolver-input.txt LPsolver-input-test1.txt

echo "Checking values of moments are as expected..."
sh moment_differ.sh $moment_filename_expected $moment_filename_test $m_thresh $u1_thresh $u2_thresh $u3_thresh $T_thresh
exit_val=$?

exit $exit_val
