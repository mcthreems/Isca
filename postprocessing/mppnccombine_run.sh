#!/usr/bin/env bash

# 1. Load the necessary tools into the environment
. /u/local/Modules/default/init/modules.sh
module purge
source $GFDL_BASE/src/extra/env/emps-hoff
module list

# Run with input file as input taken from BASH script
echo ${1}
exec ${1}/mppnccombine.x ${2}
