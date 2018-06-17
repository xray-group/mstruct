# MStruct

MStruct is a free computer program for MicroStructure analysis from powder diffraction data. 

- [https://mstruct.readthedocs.io](Documentation)
- [http://www.xray.cz/mstruct/](MSTRUCT Home Page)

## Getting Started

### Prerequisites
* [Boost](https://www.boost.org) version >=1.63
* [gsl](https://www.gnu.org/software/gsl)
* [fftw3](http://www.fftw.org)
* lapack (note: ScaLAPACK and MKL are missing a function we use)

MStruct carries the rest of the dependencies so far. This might be changed in the future.

### Compiling and Installing

Detailed installation instrucrions for different systems (Linux, Windows and macOS) can be
found on the [https://mstruct.readthedocs.io](Documentation page).

Download this repository
```bash
git clone https://github.com/xray-group/mstruct.git
```

libmstruct or mstruct can be build using the following commands:
```bash
cd mstruct/libmstruct
scons -j4 libmstruct
scons -j4 mstruct
```

libMStruct, libObjCryst and mstruct can be installed as:
```bash
scons -j4 install prefix=~/.local
```

SConscript automatically determines your system default python version and builds libMStruct
againts this python version. If you want to build libMStruct against different python version,
you can use --python-version=X.Y option for example:
```bash
scons -j4 libmstruct python-version=3.6
```
WARNING: You need Boost Python libraries to be built against the python version you want to use.

### Using

If MStruct libraries and mstruct binary were installed in non-system destinations you may like
to modify the environment in order to use MSTRUCT easily:

```bash
# let's assume an installation prefix $P
export P=~/.local
# modify environment for using sw
export PATH=$P/bin:$PATH
export LD_LIBRARY_PATH=$P/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$P/lib/python2.7/site-packages:$PYTHONPATH
# modify environment for compiling with libMStruct and libObjCryst
export CPPPATH=$P/include:$CPPPATH
export LIBRARY_PATH=$P/lib:$LIBRARY_PATH
```

Running program:
```bash
mstruct
 Beginning program ....
job type (0-data refinement,1-grid refinement)
```

Using python module:
```python
import libMstruct as mst
```
