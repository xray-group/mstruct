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
First you need to download this repository
```
git clone https://github.com/xray-group/mstruct.git
```
Then libobjcryst must be built. You can use following commands:
```
cd mstruct
cd libobjcryst/
scons -j4 lib # !modification! build->lib

```
If libobjcryst is succesfully built you can proceed building libMStruct.
```
cd ../libmstruct
scons -j4 libmstruct
```
SConscript automatically determines your system default python version and builds libMStruct againts this python version. If you want to build libMStruct against different python version, you can use --python-version=X.Y option for example:
```
scons -j4 --python-version=3.6 libmstruct
```
WARNING: You need Boost Python libraries to be built against the python version you want to use.

Next, you need to add path to your newly built libMStruct. One option is following.
*As next commands are going to modify your environment (that you may use for 
building) I suggest opening a new shell, e.g. by bash.*

```
bash
cd examples
source ./../bash-env.sh

```
You may want to check the environment

```
echo $libmstruct_path
echo $LD_LIBRARY_PATH
echo $PYTHONPATH
```
