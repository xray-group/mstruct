# MStruct
MStruct is a free computer program for MicroStructure analysis from powder diffraction data. 
## Getting Started

### Prerequisites
Boost version >=1.63
GSL
FFTW3
LAPACK
MStruct carries the rest of the dependencies so far, might be changed in the future.

### Installing

```
git clone https://github.com/xray-group/mstruct.git
cd mstruct
git checkout develop
cd libobjcryst/
scons -j4 lib # !modification! build->lib

```

```
cd ../libmstruct
scons -j4 libmstruct

```
As next commands are going to modify your environment (that you may use for 
building) I suggest opening a new shell, e.g. by bash

Open new terminal

```
bash
cd examples

```

Source the environment

```
source ./../bash-env.sh

```
You may want to check the environment

```
echo $libmstruct_path
echo $LD_LIBRARY_PATH
echo $PYTHONPATH
```
