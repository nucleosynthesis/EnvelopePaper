#!/bin/bash

#echo 'Found $SCRAM_ARCH='$SCRAM_ARCH
#echo 'Found $res='$res

if [ -z $1 ]; then
	vers=$res
	scAr=$SCRAM_ARCH
elif [ -z $2 ]; then
	vers=$1
	scAr=$SCRAM_ARCH
else
	vers=$1
	scAr=$2
fi

#echo 'Set $scAr='$scAr
#echo 'Set $vers='$vers

# SLC 6
if [ "$scAr" = "slc6_amd64_gcc472" ] || [[ $vers == *6.* ]]; then
	echo 'Setting up for SLC6'
	. /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc6-gcc47-opt/setup.sh
	. /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc6-gcc47-opt/root/bin/thisroot.sh
	export PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc47-opt/bin:$PATH
	export PYTHONPATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc47-opt/bin:$PWD:$PYTHONPATH
	export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc47-opt/lib:$PWD:$LD_LIBRARY_PATH
	export BOOSTPATH=/afs/cern.ch/sw/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc47-opt:$BOOSTPATH

# SLC 5
elif [ "$scAr" = "slc5_amd64_gcc472" ]  || [[ $vers == *5.* ]]; then
	echo 'Setting up for SLC5'
	. /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc5/setup.sh
	. /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc5-gcc47-opt/root/bin/thisroot.sh
	export PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc5-gcc47-opt/bin:$PATH
	export PYTHONPATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc5-gcc47-opt/bin:$PWD:$PYTHONPATH
	export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc5-gcc47-opt/lib:$PWD:$LD_LIBRARY_PATH
	export BOOSTPATH=/afs/cern.ch/sw/lcg/external/Boost/1.53.0_python2.7/x86_64-slc5-gcc47-opt:$BOOSTPATH

# MAC OS
elif [ "$scAr" = "osx107_amd64_gcc462" ] ; then
  echo 'Setting up for Mac OS'
  . /Applications/root/bin/thisroot.sh
  export BOOSTPATH=/Applications/boost
	export PATH=$PWD:$PATH
	export PYTHONPATH=$PWD:$PYTHONPATH
	export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
	export ROOFITSYS=$ROOTSYS

# not found
else
	# default
	echo 'Setting up for SLC6 (default)'
	. /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc6-gcc47-opt/setup.sh
	. /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc6-gcc47-opt/root/bin/thisroot.sh
	export PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc47-opt/bin:$PATH
	export PYTHONPATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc47-opt/bin:$PWD:$PYTHONPATH
	export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc47-opt/lib:$PWD:$LD_LIBRARY_PATH
	export BOOSTPATH=/afs/cern.ch/sw/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc47-opt:$BOOSTPATH
fi
