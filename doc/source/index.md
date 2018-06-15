# MStruct

## Installation

### Introduction

MStruct projects provides two main components:	
- **mstruct** program for powder diffraction data refinement
- **python module** for either interactive or automated powder diffraction analysis

Windows, MacOS and Linux are all supported. However the level of intagration is varying.
See the table below for a quick overview.


|     xxxx         | Windows       | MacOS     | Linux     |
|:--------------- |:-------------:|:---------:|:---------:|
| native binaries | available     | -         | -         |
| native source   | ask authors   | -         | available |
| Anaconda        | available     | available | available |


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

In case *an Error Message appears with the double click*
- if *libfftw3-3.dll is missing* copy/overwright *both* `mstruct.exe` and `libfftw3-3.dll` with
a pair from one of `x64` and `Win32` directories. Do not mix the pairs!
- if *MSVCP140.dll is missing* you need a set of Microsoft runtime libraries (so called *Microsoft
Visual C++ 2017 Redistributable*). It is not a compiler! You can get a version for your system
for free from [Microsoft support web](https://www.visualstudio.com/downloads/). (x86=Win32, x64=Win64)

	- [VS2017 runtime x86](https://aka.ms/vs/15/release/VC_redist.x86.exe)
	- [VS2017 runtime x64](https://aka.ms/vs/15/release/VC_redist.x64.exe)

#### Compiling with Anaconda

#### Linux native compilation

## Course - Struktura 2018

## References

