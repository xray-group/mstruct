import os
import sys
import platform

Import('env')

# Build environment configuration --------------------------------------------

# Apply CFLAGS, CXXFLAGS, LDFLAGS from the system environment.
flagnames = 'CFLAGS CPPFLAGS CXXFLAGS LDFLAGS'.split()
env.MergeFlags([os.environ.get(n, '') for n in flagnames])

# Insert LIBRARY_PATH explicitly
if env['PLATFORM'] != 'win32':
    env.PrependUnique(LIBPATH=env['ENV'].get('LIBRARY_PATH', '').split(':'))
else:
    env.PrependUnique(LIBPATH=env['ENV'].get('LIBRARY_PATH', '').split(';'))

# Insert CPPPATH explicitly
if env['PLATFORM'] != 'win32':
    env.PrependUnique(CPPPATH=env['ENV'].get('CPPPATH', '').split(':'))
else:
    env.PrependUnique(CPPPATH=env['ENV'].get('CPPPATH', '').split(';'))

# Windows specific enviroment
if env['PLATFORM'] == 'win32':
    env.PrependUnique(CPATH=env['ENV'].get('CPATH', '').split(';'))
    env.PrependUnique(CPPPATH=env['ENV'].get('CPPPATH', '').split(';'))
    env['ENV']['TMP'] = os.environ['TMP']

# Use Intel C++ compiler if requested by the user.
icpc = None
if env['tool'] == 'intelc':
    icpc = env.WhereIs('icpc')
    if not icpc:
        print("Cannot find the Intel C/C++ compiler 'icpc'.")
        Exit(1)
    env.Tool('intelc', topdir=icpc[:icpc.rfind('/bin')])

if env['PLATFORM'] != 'win32':
    fast_linkflags = ['-s']
else:
    fast_linkflags = []

# Platform specific intricacies.
if env['PLATFORM'] == 'darwin':
    env.Append(SHLINKFLAGS=['-install_name', '$TARGET.abspath'])
    env.AppendUnique(SHLINKFLAGS='-headerpad_max_install_names')
    env.AppendUnique(SHLINKFLAGS=['-undefined', 'dynamic_lookup'])
    env.AppendUnique(LINKFLAGS=['-Wl,-rpath,@loader_path/../lib'])
    fast_linkflags[:] = []

# Compiler specific options
if icpc:
    # options for Intel C++ compiler on hpc dev-intel07
    env.PrependUnique(CCFLAGS=['-w1', '-fp-model', 'precise'])
    env.PrependUnique(LIBS=['imf'])
    fast_optimflags = ['-fast', '-no-ipo']
else:
    # g++ options
    env.PrependUnique(CCFLAGS=['-Wall'])
    env.PrependUnique(CXXFLAGS=['-std=c++11'])
    fast_optimflags = ['-ffast-math']

# Configure build variants
if env['build'] == 'debug':
    env.Append(CCFLAGS='-g')
elif env['build'] == 'fast':
    env.AppendUnique(CCFLAGS=['-O3'] + fast_optimflags)
    #env.AppendUnique(CPPDEFINES={'NDEBUG' : None})
    env.AppendUnique(LINKFLAGS=fast_linkflags)

if env['profile']:
    env.AppendUnique(CCFLAGS='-pg')
    env.AppendUnique(LINKFLAGS='-pg')

# link warnings
if env['PLATFORM'] != 'win32' and env['PLATFORM'] != 'darwin':
    env.AppendUnique(LINKFLAGS=['']) # '-Wl,--' does not work with gcc4.8

# REAL=double needed for pyobjcryst
env.AppendUnique(CCFLAGS='-DREAL=double')

# required for boost with MSVC
if env['PLATFORM'] == 'win32':
    env.AppendUnique(CPPFLAGS=['/EHsc','/MD','-DBOOST_ALL_DYN_LINK','-D_DLL'])

# link flags for MSVC
if env['PLATFORM'] == 'win32':
    env.AppendUnique(LINKFLAGS=['/NODEFAULTLIB:library'])

# Lists for storing built objects and header files
env['newmatobjs'] = []
env['cctbxobjs'] = []
env['objcrystobjs'] = []
env['mstructobjs'] = []
env['binmstructobjs'] = []
env['binxmlmstructobjs'] = []
env['lib_includes'] = []
env['libmstruct_includes'] = []
env['binmstruct_includes'] = []

# Subsidiary SConscripts -----------------------------------------------------

# These will create the built objects and header file lists.
SConscript(["SConscript.cctbx", "SConscript.newmat", "SConscript.objcryst", "SConscript.mstruct"])

# Define sdist target for creating source distribution bundle
# Do so only if required to avoid extra git executions.
# Note: See .gitattributes for files excluded from sdist.
if 'sdist' in COMMAND_LINE_TARGETS:
    SConscript('SConscript.sdist')

# Top Level Targets ----------------------------------------------------------

# This retrieves the intermediate objects
newmatobjs = env["newmatobjs"]
cctbxobjs = env["cctbxobjs"]
objcrystobjs = env["objcrystobjs"]
mstructobjs = env["mstructobjs"]
binmstructobjs = env["binmstructobjs"]
binxmlmstructobjs = env["binxmlmstructobjs"]

# This builds the shared library
if env['PLATFORM'] != 'win32':
    # Additional library dependencies for ObjCryst
    ObjCrystlibs = ['lapack','m']
    ObjCrystlibpaths = []
    libobjcryst = env.SharedLibrary("libObjCryst",
            objcrystobjs + cctbxobjs + newmatobjs,
            LIBS=ObjCrystlibs) # LIBPATH=ObjCrystlibpaths
else:
    # Additional library dependencies for ObjCryst
    ObjCrystlibs = [] # Additional library dependencies for ObjCryst
    ObjCrystlibpaths = []
    env.AppendUnique(LIBPATH=[env.Dir('../../../win-prebuild/lib/x64')])
    libobjcryst = env.Library("libObjCryst",
            objcrystobjs + cctbxobjs + newmatobjs,
            LIBS=ObjCrystlibs) # LIBPATH=ObjCrystlibpaths
lib = Alias('lib', [libobjcryst, env['lib_includes']])
Default(lib)

conf = Configure(env)

# Boost configuration ----------------------------------------------------------

def find_library(possible_library_names, conf,):
    """
    Searches for library with possible different names under different versions.
    """
    proper_library_name = ''
    for name in possible_library_names:
        if conf.CheckLib(name):
            return name
    raise RuntimeError('Unable to check availability of library name (note this may also indicate other configuration problems, please check ./config.log)', possible_library_names)

boost_python_libs = []
python_libs = []

# Possible boost/numpy library names
boost_python_possible_names = ['boost_python%d%d' % env['python_version'], 'boost_python%d' % env['python_version'][0],'boost_python']
boost_numpy_possible_names = ['boost_numpy%d%d' % env['python_version'], 'boost_numpy%d' % env['python_version'][0], 'boost_numpy']
python_library_possible_names = ['python%d.%d' % env['python_version'], 'python%d.%dm' % env['python_version'], 'python%d%d' % env['python_version']]
if env['PLATFORM'] != 'darwin':
    python_libs.extend((find_library(python_library_possible_names, conf),))
boost_python_libs.extend((find_library(boost_python_possible_names, conf), find_library(boost_numpy_possible_names, conf)))

env = conf.Finish()

if env['PLATFORM'] != 'win32':
    binMStructlibs = ['fftw3', 'gsl', 'lapack', 'ObjCryst'] + ['boost_program_options',]
    libMStructlibs = ['fftw3', 'gsl', 'lapack', 'ObjCryst'] + python_libs + boost_python_libs
else:
    binMStructlibs = ['fftw3', 'gsl', 'libObjCryst', 'clapack-3.1.1-md','libf2c-3.1.1-md','blas-3.1.1-md'] + python_libs + boost_python_libs + ['boost_program_options',]
    libMStructlibs = ['fftw3', 'gsl', 'libObjCryst', 'clapack-3.1.1-md','libf2c-3.1.1-md','blas-3.1.1-md'] + python_libs + boost_python_libs + ['boost_program_options',]
MStructlibpaths = env['LIBPATH'] + ['/usr/lib', env.Dir('.')]

# This builds the shared MStruct library
libmstruct = env.SharedLibrary("libMStruct", mstructobjs, LIBS=libMStructlibs, LIBPATH=MStructlibpaths)
libms = Alias('libmstruct', [libmstruct,] + env['libmstruct_includes'])
#if env['PLATFORM'] != 'win32':
#    Depends(libmstruct, lib)

# Make sure we have @rpath/libObjCryst.dylib instead of an absolute path on Darwin
if env['PLATFORM'] == 'darwin':
    for  f in libmstruct:
        if f.get_suffix()=='.dylib':
            target_name = os.path.normpath( os.path.join('build', env['build']+'-'+env['HOST_ARCH'], f.rstr()) )
            env.AddPostAction(libmstruct,
                    'install_name_tool -change `otool -L ' + target_name + ' | grep libObjCryst.dylib | tail -n1 | tr -d "\\t" | cut -f1 -d " "` @rpath/libObjCryst.dylib ' + target_name)

# This builds mstruct binary executable
binmstruct = env.Program("mstruct_am", binmstructobjs, LIBS=binMStructlibs, LIBPATH=MStructlibpaths)
binms = Alias('mstruct', [binmstruct,] + env['binmstruct_includes'])

# This builds mstruct_xml binary executable
binxmlmstruct = env.Program("mstruct_xml", binxmlmstructobjs, LIBS=binMStructlibs, LIBPATH=MStructlibpaths)
binms_xml = Alias('mstruct_xml', [binxmlmstruct,] + env['binmstruct_includes'])

# Make sure we have @rpath/libObjCryst.dylib and @rpath/libMStruct.dylib in binaries
if env['PLATFORM'] == 'darwin':
    for build_obj, name in zip([binmstruct, binxmlmstruct],['mstruct_am','mstruct_xml']): 
        target_name = os.path.normpath( os.path.join('build', env['build']+'-'+env['HOST_ARCH'], name) )
        for libname in ['libObjCryst.dylib']:
            env.AddPostAction(build_obj,
                    'install_name_tool -change `otool -L ' + target_name + ' | grep ' + libname + ' | tail -n1 | tr -d "\\t" | cut -f1 -d " "` @rpath/' + libname + ' ' + target_name)

# Installation targets.

prefix = env['prefix']

# install-lib
libinstall = env.Install(env['libdir'], libobjcryst)
if env['PLATFORM'] == 'posix' and WhereIs('ldconfig'):
    opts = ''
    if os.getuid() != 0:  opts = '-n'
    env.AddPostAction(libinstall,
            'ldconfig %s $TARGET.dir' % opts)
if env['PLATFORM'] != 'win32':
    libinstall = env.Install(env['libdir'], libmstruct+libobjcryst)
else:
    libinstall = env.Install(env['libdir'], [f for f in libmstruct if f.get_suffix()=='.lib'])
if env['PLATFORM'] == 'posix' and WhereIs('ldconfig'):
    opts = ''
    if os.getuid() != 0:  opts = '-n'
    env.AddPostAction(libinstall,
            'ldconfig %s $TARGET.dir' % opts)
# link *.dylib -> *.so to allow python import
if env['PLATFORM'] == 'darwin':
    for  f in libmstruct:
        if f.get_suffix()=='.dylib':
            link_target = os.path.normpath( os.path.join(env['modulepath'], f.rstr()[:-5]+'so') )
            link_source = os.path.normpath( os.path.join(env['libdir'], f.rstr()) )
            env.AddPostAction(libinstall, 'rm -f ' + link_target + ' 2>/dev/null')
            env.AddPostAction(libinstall, 'ln -s ' + link_source + ' ' +  link_target)
# link site-packages/python*/*.so -> *.so to allow python import
if env['PLATFORM'] == 'posix':
    for  f in libmstruct:
        if f.get_suffix()=='.so':
            link_target = os.path.normpath( os.path.join(env['modulepath'], f.rstr()) )
            link_source = os.path.normpath( os.path.join(env['libdir'], f.rstr()) )
            env.AddPostAction(libinstall, 'rm -f ' + link_target + ' 2>/dev/null')
            env.AddPostAction(libinstall, 'ln -s ' + link_source + ' ' +  link_target)
# runtime dll for Windows
dllinstall = []
if env['PLATFORM'] == 'win32':
    dllinstall = env.Install(prefix+'/bin',  [f for f in libmstruct if f.get_suffix()=='.dll'])
    for  f in libmstruct:
        if f.get_suffix()=='.dll':
            link_source = os.path.normpath( os.path.join(env['modulepath'], f.rstr()[:-3]+'pyd') )
            link_target = os.path.normpath( os.path.join(prefix+'/bin', f.rstr()) )
            env.AddPostAction(dllinstall, 'del ' + link_source + ' 2>NUL')
            env.AddPostAction(dllinstall, 'mklink /H ' + link_source + ' ' +  link_target)
Alias('install-lib', libinstall + dllinstall)

# install-bin
if env['PLATFORM'] == 'win32':
    bininstall = env.InstallAs(prefix+'/bin/mstruct.exe', binmstruct) + env.InstallAs(prefix+'/bin/mstruct_xml.exe', binxmlmstruct)
else:
    bininstall = env.InstallAs(prefix+'/bin/mstruct', binmstruct) + env.InstallAs(prefix+'/bin/mstruct_xml', binxmlmstruct)
if env['PLATFORM'] == 'darwin':
    # DARWIN_INSTALL_NAME can be pre-set in sconscript.local
    env.SetDefault(DARWIN_INSTALL_NAME='$TARGET.abspath')
    env.AddPostAction(bininstall,
            'install_name_tool -id $DARWIN_INSTALL_NAME $TARGET')
Alias('install-bin', bininstall)

# install-python
pyinstall = []
Alias('install-modules', pyinstall)

# install-includes
ninc = len(Dir('.').path) + 1
inc_target_path = lambda f: os.path.join(env['includedir'], f.path[ninc:])
include_targets = [inc_target_path(f) for f in env['lib_includes']]
include_targets = include_targets + [inc_target_path(f) for f in env['libmstruct_includes']]
Alias('install-include', InstallAs(include_targets, env['lib_includes']+ env['libmstruct_includes']))

# install
Alias('install', ['install-lib', 'install-include', 'install-bin','install-modules'])

# vim: ft=python
