#!/bin/bash

# Grab the names of the files to be have moment values compared from
moment_filename_expected=$1
moment_filename_test=$2

# Check if thresholds have been supplied and, if so, grab their values
if [ $# -eq 2 ]; then
	m_thresh=0.000002
	u1_thresh=10^-10
	u2_thresh=10^-10
	u3_thresh=10^-10
	T_thresh=0.00003
elif [ $# -eq 7 ]; then
	m_thresh=$3
	u1_thresh=$4
	u2_thresh=$5
	u3_thresh=$6
	T_thresh=$7
else
	echo "Invalid number of arguments passed to moment_differ.sh!"
	exit 1
fi
	
# Read mass values after time-step 10
echo "- Comparing mass values..."
m_0=$(awk 'NR == 6 {print $1}' $moment_filename_expected)
m=$(awk 'NR == 6 {print $1}' $moment_filename_test)

if [ $(echo "$m - $m_0 > $m_thresh" | bc) == 1 ] || [ $(echo "$m_0 - $m > $m_thresh" | bc) == 1 ]; then
	echo "Error in mass value!"
	exit 1
fi

# Read first component of velocity values after time-step 10
echo "- Comparing first component of momentum values..."
u1_0_in=$(awk 'NR == 6 {print $2}' $moment_filename_expected)
u1_0=$(echo "$u1_0_in" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
u1_in=$(awk 'NR == 6 {print $2}' $moment_filename_test)
u1=$(echo "$u1_in" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')

if [ $(echo "$u1 - $u1_0 > $u1_thresh" | bc -l) == 1 ] || [ $(echo "$u1_0 - $u1 > $u1_thresh" | bc -l) == 1 ]; then
	echo "Error in first component of momentum!"
	exit 1
fi

# Read second component of velocity values after time-step 10
echo "- Comparing second component of momentum values..."
u2_0_in=$(awk 'NR == 6 {print $3}' $moment_filename_expected)
u2_0=$(echo "$u2_0_in" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
u2_in=$(awk 'NR == 6 {print $3}' $moment_filename_test)
u2=$(echo "$u2_in" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')

if [ $(echo "$u2 - $u2_0 > $u2_thresh" | bc -l) == 1 ] || [ $(echo "$u2_0 - $u2 > $u2_thresh" | bc -l) == 1 ]; then
	echo "Error in second component of momentum!"
	exit 1
fi

# Read third component of velocity values after time-step 10
echo "- Comparing third component of momentum values..."
u3_0_in=$(awk 'NR == 6 {print $4}' $moment_filename_expected)
u3_0=$(echo "$u3_0_in" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
u3_in=$(awk 'NR == 6 {print $4}' $moment_filename_test)
u3=$(echo "$u3_in" | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')

if [ $(echo "$u3 - $u3_0 > $u3_thresh" | bc -l) == 1 ] || [ $(echo "$u3_0 - $u3 > $u3_thresh" | bc -l) == 1 ]; then
	echo "Error in third component of momentum!"
	exit 1
fi

# Read total energy values after time-step 10
echo -e "- Comparing total energy values... \n"
T_0=$(awk 'NR == 6 {print $9}' $moment_filename_expected)
T=$(awk 'NR == 6 {print $9}' $moment_filename_test)
if [ $(echo "$T - $T_0 > $T_thresh" | bc) == 1 ] || [ $(echo "$T_0 - $T > $u2_thresh" | bc) == 1 ]; then
	echo "Error in total energy!"
	exit 1
fi

exit 0
