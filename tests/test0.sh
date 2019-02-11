#!/bin/bash

mv LPsolver-input-Test0.txt LPsolver-input.txt
../source/solver
mv LPsolver-input.txt LPsolver-input-Test0.txt

sh moment_differ.sh
exit_val=$?

exit $exit_val
