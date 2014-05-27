#!/bin/bash

version=`cat /proc/version | awk -F " " '{print $3}'`
if version=='2.6.32-431.17.1.el6.x86_64'; then
	# SLC 6
	echo 'Setting up for SLC6'
	. /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc6-gcc47-opt/setup.sh
	. /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc6-gcc47-opt/root/bin/thisroot.sh
	export PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc47-opt/bin:$PATH
	export PYTHONPATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc47-opt/bin:$PWD:$PYTHONPATH
	export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc47-opt/lib:$PWD:$LD_LIBRARY_PATH
	export BOOSTPATH=/afs/cern.ch/sw/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc47-opt:$BOOSTPATH
else
	# SLC 5
	echo 'Setting up for SLC5'
	. /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc5/setup.sh
	. /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc5-gcc47-opt/root/bin/thisroot.sh
	export PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc5-gcc47-opt/bin:$PATH
	export PYTHONPATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc5-gcc47-opt/bin:$PWD:$PYTHONPATH
	export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc5-gcc47-opt/lib:$PWD:$LD_LIBRARY_PATH
	export BOOSTPATH=/afs/cern.ch/sw/lcg/external/Boost/1.53.0_python2.7/x86_64-slc5-gcc47-opt:$BOOSTPATH
fi
