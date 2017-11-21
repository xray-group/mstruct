#!/bin/bash

# DO NOT MOVE THIS SCRIPT
# source it as: source /path/to/the/script/script.name

SDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export libobjcryst_path=$SDIR/../libobjcryst/build/fast-x86_64/
export libmstruct_path=$SDIR/../libmstruct/build/fast-x86_64/

export LD_LIBRARY_PATH=$libobjcryst_path:$libmstruct_path:$LD_LIBRARY_PATH
export CPATH=$libobjcryst_path:$libobjcryst_path/ObjCryst:$libmstruct_path:$libmstruct_path/MStruct:$CPATH:
export LIBRARY_PATH=$libobjcryst_path:$libmstruct_path:$LIBRARY_PATH
