# MStruct

software/library for MicroStructure analysis by powder diffraction 

## Installation

### MStruct and MStruct GUI

MStruct software for micro-structure analysis from powder diffraction data has two parts nowadyas:
- **MStruct-cli** - the original command line interface extended with Python API
- **MStruct GUI** - a Java based user intrface using the **MStruct-cli** as a computing backend

Credit for **MStruct GUI** goes to Jakub Vojtisek, Milan Dopita and Lukas Horak (MFF, Charles University, Prague).

The binary distributions for the MStruct software can be downloaded from the
[shared Dropbox folder](https://www.dropbox.com/scl/fo/4v3m7v401vg689frq0lnw/AMe6l_ygImbikjbyuymWUyg?rlkey=t7qjttvlrwd51ij5t0iclx3ck&dl=0) (MStructGUI-0.4 + mstruct-0.15.10).

| OS                    | name                            | md5                              |
|:----------------------|:-------------------------------:|:--------------------------------:|
| Windows (x64)         | [MStruct Installer-0.4 (x64).msi](https://www.dropbox.com/scl/fi/93nbz9v9qfg8ncd16stle/MStruct-Installer-0.4-x64.msi?rlkey=7g4v4a65rcbp8931pimljmmzv&dl=0) | 1e2950e959dc83f8717c60700c23a28a |
| Linux (ubuntu, amd64) | [mstruct-gui_0.4_amd64.deb](https://www.dropbox.com/scl/fi/uwm8e9kl70gdznzsis4l6/mstruct-gui_0.4_amd64.deb?rlkey=ty7uv53h46ietwn9qvn51k7yw&dl=0)       | -                                |
| Linux (rhel, amd64)   | [mstruct-gui-0.4-1.x86_64.rpm](https://www.dropbox.com/scl/fi/l271ptfo01xgtdef1lnfr/mstruct-gui-0.4-1.x86_64.rpm?rlkey=oka933c5w7mzr6jazzqf5jtpc&dl=0)    | -                                |
| macOS (amd64)         | [MStruct GUI - 0.4(amd64).pkg](https://www.dropbox.com/scl/fi/u7912fykj40rq5pkrn4kt/MStruct-GUI-0.4-amd64.pkg?rlkey=bpd86h3849m2gjv2jp3w6iht6&dl=0)    | -                                |
| macOS (aarch64)       | [MStruct GUI - 0.4(aarch64).pkg](https://www.dropbox.com/scl/fi/8ngabqimf4sodc596zu44/MStruct-GUI-0.4-aarch64.pkg?rlkey=je9g3eqdw5j078wyhh5523b0f&dl=0)  | -                                |

The first installation step is running the intaller. It is strongly recomended to install into a default location.

Installers provide both the CLI and GUI. One starter is named *GUI*. The other one like *Prompt*. The *Prompt* provides a command line terminal with MStruct CLI and Python module available.

MStruct GUI requires for its function a configuration file **MStructGUI.properties**. This file contains beside others the path to the MStruct CLI. It may be good to copy this file with some additional resources into the *user space*. Please follow the OS specific instructions below for this step.

#### Windows

Windows installer copies the required resources in Documents folder and the starter icon on the Desktop starts MStruct GUI in that folder. So no additional actions are required.

##### Know issues

MStruct GUI does not function properly on localized (non-English) Windows. However it is quite easy to switch between languages in Win 10/11. Type "Language set..." in the search field to get to the Windows language settings. A new login will be requred after changing the language to English.

An alternative is to use Windows Subsystem for Linux (WSL) and install one of the Linux packeges.

#### Linux

The default Linux installation location is `/opt/xray-group/mstruct`. The user specific resources can be found in the `shared` folder. So ideally copy them to a user writable space and start MStruct GUI.

```bash
cp -R /opt/xray-group/mstruct/shared/mstruct-gui .   # copy resources
cd msttruct-gui                                      # change working directory
mstruct-gui                                          # start mstruct gui
````

Beside GUI the installation provides an environment script with mstruct CLI and Python API available:

```bash
source /opt/xray-group/mstruct/source_me.bash
mstruct_xml --version
````

#### macOS

The default macOS installation location is `/Applications/MStruct_GUI`. The user specific resources can be found in the `shared` folder. So ideally copy them to a user writable space and start MStruct GUI.

```bash
cp -R /Applications/MStruct_GUI/shared/mstruct-gui .  # copy resources
cd msttruct-gui                                       # change working directory
mstruct-gui                                           # start mstruct gui
````

Beside GUI the installation provides an environment script with mstruct CLI and Python API available:

```bash
source /Applications/MStruct_GUI/source_me.bash
mstruct_xml --version
````

The MStruct GUI or an mstruct terminal can be also started by double clicking on the command-starters. However remmeber to copy **MStructGUI.properties** to your home folder for a proper function of MStruct GUI.

##### Know issues

Rescaling MStruct GUI window with modal dialogues may make it cumberstone to work with MStruct. It may help to change the *Dock* System Settings.

On newer (Sequoia v15) macOS versions go to *System Settings* and type *Desktop & Dock*. Scroll down to *Windows* section and select **Prefer tabs when opening documents** **Never**. On more earlier versions (Monterey v12) look for this option in the *General* Tab.

### MStruct CLI introduction

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
- **Windows binaries** are small (few MBs), easy to download and run.
- **Anaconda** is the most universal. You will get all MStruct features without
  limitations but you need to have Anaconda enviroment. Anaconda occupies around
  2-3 GBs. However you may use it also for something else. You will need to
  compile MStruct yourself, short instructions are provided, we tested it but
  some issues are hard to exluclude with any effort (just try it).
  There is a less demanding **Miniconda** alternative.
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
conda install python=3.11 boost=1.78 lapack fftw gsl scons bzip2 git blas=*=*mkl gcc gxx

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

## Course - Struktura 2025

### Instructions

You should have `mstruct` binary running on your laptop. It should be pretty straighforward
for Windows, unfortunately for macOS or Linux you need to compile the source. *No stress in
case you fail, we will make it working at the beginning of the course.*

### (prelliminary) Time plan

- Mon Sep 8th, 2025, 17:50 - 18:30. **Installation of MStruct GUI** (Zdenek, Milan, Lukas)
- Mon Sep 8th, 2025, 18:30 - 19:30. **MStruct overview: The basics** (optional talk by Zdenek)
- Mon Sep 8th, 2025, 20:30 - 21:30. **Real structure analysis: The basics** (Milan, Lukas Zdenek)
	  - *Residual stress in thin TiO2 films*: residual stress and refraction correction
    - *Nanocrystallie TiO2 powders for catalytic applications*: spherical crystallites,
	crystallite size distribution, isotropic phenomenological microstrain, crystal
	structure parameters, quantitative phase analysis
    - *Instrumental profile - "Fitovacka"* (Lukas)
    - *Python API* (Zdenek)
- Mon Sep 8th, 2025, 21:30 - 19:40. **Challenging samples: Advanced analysis** (Milan)
    - *Ultrathin Pt nanocrystalline films*: refraction correction, thin film correction,
	microstrain, stacking faults
	  - *Copper-Gold nano-spheres*: Using individual peak parameters with instrumental
	correction
- Fri Sep 12, 2025, 9:30 -- (optional): **Hands on data** (Milan, Lukas, Zdenek)

## References

1. [MSTRUCT materials and tutorial data](https://github.com/xray-group/mstruct-materials)
2. [MSTRUCT Home page](https://www.xray.cz/mstruct/)
3. [Basic MSTRUCT Tutorial](https://www.xray.cz/mstruct/mstruct-basic-ex.html)
4. [MSTRUCT prezentation (2013)](https://www.xray.cz/mstruct/mstruct-pres-2013.pdf)
5. [MSTRUCT on GitHub](https://github.com/xray-group/mstruct)
6. [DOI: 10.1017/S0885715614000852](http://dx.doi.org/10.1017/S0885715614000852)