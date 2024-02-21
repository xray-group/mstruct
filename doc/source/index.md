# MStruct

software/library for MicroStructure analysis by powder diffraction 

## Installation

### Introduction

MStruct projects provides two main components:	
- **mstruct** program for powder diffraction data refinement
- **python module** for either interactive or automated powder diffraction analysis

Windows, macOS and Linux are all supported. However the level of intagration is varying.
See the table below for a quick overview.

|                 | Windows       | macOS     | Linux     |
|:--------------- |:-------------:|:---------:|:---------:|
| native binaries | available     | -         | available |
| native source   | depricated    | -         | available |
| Anaconda        | available     | available | available |

<!--
|                 | Windows       | macOS     | Linux     |
| native binaries | available     | -         | -         |
| native source   | ask authors   | -         | available |
| Anaconda        | available     | available | available |
!-->

In short:
- **Windows binaries** are small (few MBs), easy to download and run but you will
  miss the Python module. Maybe you do not care.
- **Anaconda** is the most universal. You will get all MStruct features without
  limitations but you need to have Anaconda enviroment. Anaconda occupies around
  2-3 GBs. However you may use it also for something else. You will need to
  compile MStruct yourself, short instructions are provided, we tested it but
  some issues are hard to exluclude with any effort (just try it).
- **Linux native compilation** will give you all features and will not use
  much space. Compilation process is similar to Anaconda. Use of Anaconda is still
  adviced mainly in order to protect your system against dirty MStruct features :-)

### Installation instructions

#### Windows binaries

1. *Download* the latest version of [mstruct-win.zip](http://www.xray.cz/mstruct/mstruct-win.zip)
   from the [Czecho-Slovak Crystallographic Server](http://www.xray.cz/)
2. *Unzip*
3. *Open* `mstruct-win` folder, where you can find
	- the program `mstruct.exe` and some configuration and example files
	- subdirectories `x64` and `Win32` with the program for 32 or 64 bit Windows
4. *Double click* `mstruct.exe` if you get the console windows with some text the software
   is working

You can copy or move `mstruct.exe` as you like but remeber to *keep libfftw3-3.dll together with
mstruct.exe* in the same directory!

**Troubleshouting**

If *an Error Message appears with the double click*
- in case *libfftw3-3.dll is missing* copy/overwright *both* `mstruct.exe` and `libfftw3-3.dll` with
an appropriate pair from one of `x64` and `Win32` directories. Do not mix the pairs!
- in case *MSVCP140.dll is missing* you need a set of Microsoft runtime libraries (so called *Microsoft
Visual C++ 2017 Redistributable*). It is not a compiler! You can get a version for your system
for free from [Microsoft support web](https://www.visualstudio.com/downloads/). (x86=Win32, x64=Win64)

	- [VS2017 runtime x86](https://aka.ms/vs/15/release/VC_redist.x86.exe)
	- [VS2017 runtime x64](https://aka.ms/vs/15/release/VC_redist.x64.exe)

#### Obtaining source

MStruct is a free open source software. The source code can be obtained from the public
[MSTRUCT GitHub repo](https://github.com/xray-group/mstruct).

```bash
# clone it
git clone https://github.com/xray-group/mstruct.git

# or download as zip, manually: top right green button
```

#### Compiling with Anaconda

##### Installing Anaconda

[Anaconda](https://www.anaconda.com) is a popular Python data science platform and
scientific software for personal computers with Windows, macOS or Linux.

The most straighforward way is to get a graphical [installer](https://www.anaconda.com/download/)
- no need to sign (web download)
- administrative rights not required for a personal installation
- no need to add to PATH (installation option)
- preferred choice to Register Anaconda as a system Python (installation option), good
choice if you do not have any or do not care


**Anaconda Navigator** in Windows menu:

![Anaconda Navigator in Windows menu](https://raw.githubusercontent.com/xray-group/mstruct/doc/doc/source/figs/installation-anaconda-menu.png "Anaconda Navigator in Windows menu")

**Creating Anaconda environment** with name='mst':

![Creating Anaconda environment](https://raw.githubusercontent.com/xray-group/mstruct/doc/doc/source/figs/installation-anaconda-create-env.png "Creating Anaconda environment")

**Activating Anaconda environment** name='mst':

![Activating Anaconda environment](https://raw.githubusercontent.com/xray-group/mstruct/doc/doc/source/figs/installation-anaconda-open-env.png "Activating Anaconda environment")

**Installing sw in Conda environment**:

![Installing sw with Anaconda](https://raw.githubusercontent.com/xray-group/mstruct/doc/doc/source/figs/installation-anaconda-conda-install.png "Installing sw with Anaconda")

##### Windows with Anaconda

For Windows **Python3 (x64) is strongly adviced!**

```bash
# add 'conda-forge' channel
conda config --add channels conda-forge
# install required packages
conda install python=3.11 boost=1.78 lapack fftw gsl scons bzip2 git blas=*=*mkl

# git clone or download ZIP
# git clone https://github.com/xray-group/mstruct.git
# wget https://github.com/xray-group/mstruct/archive/r0.15.zip

# swith to project directory
cd mstruct/libmstruct

# set prefix path %P% where your Anaconda environment is installed
# do not forget the name='mst' at the end
set P=C:/..../Anaconda3/envs/mst
# alternatively
set P=conda info --json | grep default_prefix | cut -d '"' -f4
# set environment variables
set CPPPATH=%P%/Library/include;%P%/include
set LIBRARY_PATH=%P%/libs;%P%/Library/lib

# build and install everything
scons -j4 install prefix=%P%/Library modulepath=%P%/Lib/site-packages
```

**Test**

```bash
# type "mstruct"
mstruct
# you should see text (CTRL+C to exit)
 Beginning program ....
job type (0-data refinement,1-grid refinement)
```

##### macOS with Conda

```bash
# setup environment, e.g. name='mst'
# see instructions using Anaconda Navigator for Windows
# or cmd-line instruction for Linux

# activate environment name='mst'
source activate mst

# install required packages
conda install python=3.11 boost=1.78 lapack fftw gsl scons bzip2

# resolve the prefix directory P of the active Anaconda environment
P="$(conda info --json | grep default_prefix | cut -d\" -f4)"
# use it to setup environment variables
export CPPPATH=$P/include
export LIBRARY_PATH=$P/lib
export LDFLAGS=-Wl,-rpath,$P/lib

# git clone or download ZIP
# git clone https://github.com/xray-group/mstruct.git
# wget https://github.com/xray-group/mstruct/archive/r0.15.zip

# swith to project directory
cd mstruct/libmstruct

# build library
scons -j4 libmstruct

# build mstruct
scons -j4 mstruct

# (optional) build and install everything   
scons -j4 install prefix=$P
```

**Test**

```bash
# type "mstruct"
mstruct
# you should see text (CTRL+C to exit)
 Beginning program ....
job type (0-data refinement,1-grid refinement)
```

Note unfortunatelly you need to set the `LD_LIBRARY_PATH` every time you activate
the environment.

##### Linux with Conda

```bash
# (optional block)
# ----------------
# see installed environments
conda info --envs
# setting up an evironment name='mst'
conda create -n mst
# activating evironment name='mst'
source activate mst

# gls, fftw3, lapack, scons are required
conda install -c conda-forge python=3.11 boost=1.78 lapack fftw gsl scons bzip2 blas=*=*mkl
# note:
# - check if you like the python version and if it is consitent across the packages
# - prefer builds that are close to your current environment

# get source, git-clone or download and unpack the source
git clone https://github.com/xray-group/mstruct.git

# swith to project directory
cd mstruct/libmstruct

# resolve the prefix directory P of the active Anaconda environment
P="$(conda info --json | grep default_prefix | cut -d\" -f4)"
export CPPPATH=$P/include
export LIBRARY_PATH=$P/lib
export LD_LIBRARY_PATH=$P/lib

# build library
scons -j4 libmstruct

# build mstruct
scons -j4 mstruct

# (optional) build and install everything   
scons -j4 install prefix=$P
```

**Test**

```bash
# type "mstruct"
mstruct
# you should see text (CTRL+C to exit)
 Beginning program ....
job type (0-data refinement,1-grid refinement)
```

```bash
P="$(conda info --json | grep default_prefix | cut -d\" -f4)"
export LD_LIBRARY_PATH=$P/lib
```

##### Conda tips and tricks

tricks from Honza

#### Linux native compilation

```bash
# (optional) gls, fftw3, lapack, python and scons are required
sudo apt-get install libgsl-dev libfftw3-dev liblapack-dev libpython3.11-dev scons

# (optional) boost>=1.63 is required, we may want to use a specific one
export B=~/sw/boost_1_78_0
export CPPPATH=$B/include:$CPPPATH
export LIBRARY_PATH=$B/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=$B/lib:$LD_LIBRARY_PATH

# (optional) prepare user env for installation (define prefix)
export P=~/.local
# (optional) make sure we have place for python modules
mkdir -p $P/lib/python3.11/site-packages

# build library
scons -j4 libmstruct

# build mstruct
scons -j4 mstruct

# (optional) build and install everything
scons -j4 install prefix=$P

# (optional) we may want to activate the installation
export PATH=$P/bin:$PATH
export PYTHONPATH=$P/lib/python3.11/site-packages:$PYTHONATH
export LD_LIBRARY_PATH=$P/lib:$LD_LIBRARY_PATH
```

**Test**

```bash
# type "mstruct" (if you activated the installation as indicated above)
mstruct
# you should see text (CTRL+C to exit)
 Beginning program ....
job type (0-data refinement,1-grid refinement)
```

#### Text editors

some advices to text editors

#### Plotting tools

links to plotting tools

## Course - Struktura 2018

### Instructions

You should have `mstruct` binary running on your laptop. It should be pretty straighforward
for Windows, unfortunately for macOS or Linux you need to compile the source. *No stress in
case you fail, we will make it working at the beginning of the course.*

### Time plan

- Mon June 18, 2018, 20:00 - 21:30. **Real structure analysis: The basics** (Zdenek)
    - *Nanocrystallie TiO2 powders for catalytic applications*: spherical crystallites,
	crystallite size distribution, isotropic phenomenological microstrain, crystal
	structure parameters, quantitative phase analysis
	- *Residual stress in thin TiO2 films*: residual stress and refraction correction

- Tue June 19, 2018, 17:10 - 19:40. **Challenging samples: Advanced analysis** (Milan)
    - *Ultrathin Pt nanocrystalline films*: refraction correction, thin film correction,
	microstrain, stacking faults
	- *Copper-Gold nano-spheres*: Using individual peak parameters with instrumental
	correction

- Wed June 20, 2018, 17:30 - 19:30. (optional): **Hands on own data** (Milan, Zdenek)

- Thu June 20, 2018, 14:00 --. (optional): **Individual discussions** (Zdenek)
	- *amorphous content determination*
	- *texture*

## References

1. [MSTRUCT Home page](http://www.xray.cz/mstruct/)
2. [Basic MSTRUCT Tutorial](http://www.xray.cz/mstruct/mstruct-basic-ex.html)
3. [MSTRUCT prezentation (2013)](http://www.xray.cz/mstruct/mstruct-pres-2013.pdf)
4. [MSTRUCT on GitHub](https://github.com/xray-group/mstruct)
5. [DOI: 10.1017/S0885715614000852](http://dx.doi.org/10.1017/S0885715614000852)



