#!/usr/bin/env python

# **************************************************************************
# *
# * Authors:     I. Foche Perez (ifoche@cnb.csic.es)
# *              J. Burguet Castell (jburguet@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'ifoche@cnb.csic.es'
# *
# **************************************************************************

# Builders and pseudobuilders used be SConscript to install things.


import os
from os.path import join, abspath
import platform
import SCons.Script


# OS boolean vars
MACOSX = (platform.system() == 'Darwin')
WINDOWS = (platform.system() == 'Windows')
LINUX = (platform.system() == 'Linux')

# URL where we have most of our tgz files for libraries, modules and packages.
URL_BASE = 'http://scipionwiki.cnb.csic.es/files/scipion/software'

# Define our builders.
Download = Builder(action='wget -nv $SOURCE -c -O $TARGET')
Untar = Builder(action='tar -C $cdir --recursive-unlink -xzf $SOURCE')

# Create the environment the whole build will use.
env = Environment(ENV=os.environ,
                  tools=['Make', 'AutoConfig'],
                  toolpath=[join('software', 'install', 'scons-tools')])

# Message from autoconf and make, so we don't see all its verbosity.
env['AUTOCONFIGCOMSTR'] = "Configuring $TARGET from $SOURCES"
env['MAKECOMSTR'] = "Compiling & installing $TARGET from $SOURCES "

# Add the path to dynamic libraries so the linker can find them.
if LINUX:
    env.AppendUnique(LIBPATH=os.environ.get('LD_LIBRARY_PATH'))
elif MACOSX:
    print "OS not tested yet"
    env.AppendUnique(LIBPATH=os.environ.get('DYLD_FALLBACK_LIBRARY_PATH'))
elif WINDOWS:
    print "OS not tested yet"


#  ************************************************************************
#  *                                                                      *
#  *                           Pseudobuilders                             *
#  *                                                                      *
#  ************************************************************************

# We have 4 "Pseudo-Builders" http://www.scons.org/doc/HTML/scons-user/ch20.html
#
# They are:
#   addLibrary    - install a library
#   addModule     - install a Python module
#   addPackage    - install an EM package
#   manualInstall - install by manually running commands
#
# Their structure is similar:
#   * Define reasonable defaults
#   * Set --with-<name> option as appropriate
#   * Concatenate builders
#
# For the last step we concatenate the builders this way:
#   target1 = Builder1(env, target1, source1)
#   SideEffect('dummy', target1)
#   target2 = Builder2(env, target2, source2=target1)
#   SideEffect('dummy', target2)
#   ...
#
# So each target becomes the source of the next builder, and the
# dependency is solved. Also, we use SideEffect('dummy', ...) to
# ensure it works in a parallel build (see
# http://www.scons.org/wiki/SConsMethods/SideEffect), and does not try
# to do one step while the previous one is still running in the background.

def addLibrary(env, name, tar=None, buildDir=None, targets=None,
               url=None, flags=[], autoConfigTarget='Makefile',
               deps=[], clean=[], default=True):
    """Add library <name> to the construction process.

    This pseudobuilder downloads the given url, untars the resulting
    tar file, configures the library with the given flags, compiles it
    (in the given buildDir) and installs it. It also tells SCons about
    the proper dependencies (deps).

    If default=False, the library will not be built unless the option
    --with-<name> is used.

    Returns the final targets, the ones that Make will create.

    """
    # Use reasonable defaults.
    tar = tar or ('%s.tgz' % name)
    url = url or ('%s/external/%s' % (URL_BASE, tar))
    buildDir = buildDir or tar.rsplit('.tar.gz', 1)[0].rsplit('.tgz', 1)[0]
    targets = targets or ['lib/lib%s.so' % name]

    # Add "software/lib" and "software/bin" to LD_LIBRARY_PATH and PATH.
    def pathAppend(var, value):
        valueOld = os.environ.get(var, '')
        i = 0  # so if flags is empty, we put it at the beginning too
        for i in range(len(flags)):
            if flags[i].startswith('%s=' % var):
                valueOld = flags.pop(i).split('=', 1)[1] + ':' + valueOld
                break
        flags.insert(i, '%s=%s:%s' % (var, value, valueOld))
    pathAppend('LD_LIBRARY_PATH', abspath('software/lib'))
    pathAppend('PATH', abspath('software/bin'))

    # Install everything in the appropriate place.
    flags += ['--prefix=%s' % abspath('software'),
              '--libdir=%s' % abspath('software/lib')]  # not lib64

    # Add the option --with-name, so the user can call SCons with this
    # to activate the library even if it is not on by default.
    if not default:
        AddOption('--with-%s' % name, dest=name, action='store_true',
                  help='Activate library %s' % name)

    # Create and concatenate the builders.
    tDownload = Download(env, 'software/tmp/%s' % tar, Value(url))
    SideEffect('dummy', tDownload)  # so it works fine in parallel builds
    tUntar = Untar(env, 'software/tmp/%s/configure' % buildDir, tDownload,
                   cdir='software/tmp')
    SideEffect('dummy', tUntar)  # so it works fine in parallel builds
    Clean(tUntar, 'software/tmp/%s' % buildDir)
    tConfig = env.AutoConfig(
        source=Dir('software/tmp/%s' % buildDir),
        AutoConfigTarget=autoConfigTarget,
        AutoConfigSource='configure',
        AutoConfigParams=flags,
        AutoConfigStdOut='software/log/%s_config.log' % name)
    SideEffect('dummy', tConfig)  # so it works fine in parallel builds
    tMake = env.Make(
        source=tConfig,
        target=['software/%s' % t for t in targets],
        MakePath='software/tmp/%s' % buildDir,
        MakeEnv=os.environ,
        MakeTargets='install',
        MakeStdOut='software/log/%s_make.log' % name)
    SideEffect('dummy', tMake)  # so it works fine in parallel builds

    # Clean the special generated files
    for cFile in clean:
        Clean(tMake, cFile)

    # Add the dependencies.
    for dep in deps:
        Depends(tConfig, dep)

    if default or GetOption(name):
        Default(tMake)

    return tMake


def addModule(env, name, tar=None, buildDir=None, targets=None,
              url=None, flags=[], deps=[], clean=[], default=True):
    """Add Python module <name> to the construction process.

    This pseudobuilder downloads the given url, untars the resulting
    tar file, configures the module with the given flags, compiles it
    (in the given buildDir) and installs it. It also tells SCons about
    the proper dependencies (deps).

    If default=False, the module will not be built unless the option
    --with-<name> is used.

    Returns the final target (software/lib/python2.7/site-packages/<name>).

    """
    # Use reasonable defaults.
    tar = tar or ('%s.tgz' % name)
    url = url or ('%s/python/%s' % (URL_BASE, tar))
    buildDir = buildDir or tar.rsplit('.tar.gz', 1)[0].rsplit('.tgz', 1)[0]
    targets = targets or [name]
    flags += ['--prefix=%s' % abspath('software')]

    # Add the option --with-name, so the user can call SCons with this
    # to activate the module even if it is not on by default.
    if not default:
        AddOption('--with-%s' % name, dest=name, action='store_true',
                  help='Activate module %s' % name)

    # Create and concatenate the builders.
    tDownload = Download(env, 'software/tmp/%s' % tar, Value(url))
    SideEffect('dummy', tDownload)  # so it works fine in parallel builds
    tUntar = Untar(env, 'software/tmp/%s/setup.py' % buildDir, tDownload,
                   cdir='software/tmp')
    SideEffect('dummy', tUntar)  # so it works fine in parallel builds
    Clean(tUntar, 'software/tmp/%s' % buildDir)
    tInstall = env.Command(
        ['software/lib/python2.7/site-packages/%s' % t for t in targets],
        tUntar,
        Action('PYTHONHOME="%(root)s" LD_LIBRARY_PATH="%(root)s/lib" '
               'PATH="%(root)s/bin:%(PATH)s" '
               '%(root)s/bin/python setup.py install %(flags)s > '
               '%(root)s/log/%(name)s.log 2>&1' % {'root': abspath('software'),
                                                   'PATH': os.environ['PATH'],
                                                   'flags': ' '.join(flags),
                                                   'name': name},
               'Compiling & installing %s > software/log/%s.log' % (name, name),
               chdir='software/tmp/%s' % buildDir))
    SideEffect('dummy', tInstall)  # so it works fine in parallel builds

    # Clean the special generated files
    for cFile in clean:
        Clean(lastTarget, cFile)

    # Add the dependencies.
    for dep in deps:
        Depends(tInstall, dep)

    if default or GetOption(name):
        Default(tInstall)

    return tInstall


def addPackage(env, name, tar=None, buildDir=None, url=None,
               extraActions=[], deps=[], clean=[], default=True):
    """Add external (EM) package <name> to the construction process.

    This pseudobuilder downloads the given url, untars the resulting
    tar file and copies its content from buildDir into the
    installation directory <name>. It also tells SCons about the
    proper dependencies (deps).

    extraActions is a list of (target, command) that should be
    executed after the package is properly installed.

    If default=False, the package will not be built unless the option
    --with-<name> or --with-all-packages is used.

    Returns the final target (software/em/<name>).

    """
    # Use reasonable defaults.
    tar = tar or ('%s.tgz' % name)
    url = url or ('%s/em/%s' % (URL_BASE, tar))
    buildDir = buildDir or tar.rsplit('.tar.gz', 1)[0].rsplit('.tgz', 1)[0]

    # Add the option --with-<name>, so the user can call SCons with this
    # to get the package even if it is not on by default.
    AddOption('--with-%s' % name, dest=name, metavar='%s_HOME' % name.upper(),
              nargs='?', const='unset',
              help=("Get package %s. With no argument, download and "
                    "install it. To use an existing installation, pass "
                    "the package's directory." % name))
    # So GetOption(name) will be...
    #   None      if we did *not* pass --with-<name>
    #   'unset'   if we passed --with-<name> (nargs=0)
    #   PKG_HOME  if we passed --with-<name>=PKG_HOME (nargs=1)

    # See if we have used the --with-<package> option and exit if appropriate.
    if GetOption('withAllPackages'):
        defaultPackageHome = 'unset'
        # we asked to install all packages, so it is at least as if we
        # also did --with-<name>
    else:
        defaultPackageHome = None
        # by default it is as if we did not use --with-<name>

    packageHome = GetOption(name) or defaultPackageHome

    # If we do have a local installation, link to it and exit.
    if packageHome != 'unset':  # default value when calling only --with-package
        # Just link to it and do nothing more.
        return env.Command(
            Dir('software/em/%s/bin' % name),
            Dir(packageHome),
            Action('rm -rf %s && ln -v -s %s %s' % (name, packageHome, name),
                   'Linking package %s to software/em/%s' % (name, name),
                   chdir='software/em'))

    # Donload, untar, link to it and execute any extra actions.
    tDownload = Download(env, 'software/tmp/%s' % tar, Value(url))
    SideEffect('dummy', tDownload)  # so it works fine in parallel builds
    tUntar = Untar(env, Dir('software/em/%s/bin' % buildDir), tDownload,
                   cdir='software/em')
    SideEffect('dummy', tUntar)  # so it works fine in parallel builds
    Clean(tUntar, 'software/em/%s' % buildDir)
    if buildDir != name:
        # Yep, some packages untar to the same directory as the package
        # name (hello Xmipp), and that is not so great. No link to it.
        tLink = env.Command(
            'software/em/%s/bin' % name,  # TODO: find smtg better than "/bin"
            Dir('software/em/%s' % buildDir),
            Action('rm -rf %s && ln -v -s %s %s' % (name, buildDir, name),
                   'Linking package %s to software/em/%s' % (name, name),
                   chdir='software/em'))
    else:
        tLink = tUntar  # just so the targets are properly connected later on
    SideEffect('dummy', tLink)  # so it works fine in parallel builds
    lastTarget = tLink

    for target, command in extraActions:
        lastTarget = env.Command('software/em/%s/%s' % (name, target),
                                 lastTarget,
                                 Action(command, chdir='software/em/%s' % name))
        SideEffect('dummy', lastTarget)  # so it works fine in parallel builds

    # Clean the special generated files
    for cFile in clean:
        Clean(lastTarget, cFile)
    # Add the dependencies. Do it to the "link target" (tLink), so any
    # extra actions (like setup scripts) have everything in place.
    for dep in deps:
        Depends(tLink, dep)

    if default or packageHome:
        Default(lastTarget)

    return lastTarget


def manualInstall(env, name, tar=None, buildDir=None, url=None,
                  extraActions=[], deps=[], clean=[], default=True):
    """Just download and run extraActions.

    This pseudobuilder downloads the given url, untars the resulting
    tar file and runs extraActions on it. It also tells SCons about
    the proper dependencies (deps).

    extraActions is a list of (target, command) that should be
    executed after the package is properly installed.

    If default=False, the package will not be built unless the option
    --with-<name> is used.

    Returns the final target in extraActions.

    """
    # Use reasonable defaults.
    tar = tar or ('%s.tgz' % name)
    url = url or ('%s/external/%s' % (URL_BASE, tar))
    buildDir = buildDir or tar.rsplit('.tar.gz', 1)[0].rsplit('.tgz', 1)[0]

    # Add the option --with-name, so the user can call SCons with this
    # to activate it even if it is not on by default.
    if not default:
        AddOption('--with-%s' % name, dest=name, action='store_true',
                  help='Activate %s' % name)

    # Donload, untar, and execute any extra actions.
    tDownload = Download(env, 'software/tmp/%s' % tar, Value(url))
    SideEffect('dummy', tDownload)  # so it works fine in parallel builds
    tUntar = Untar(env, 'software/tmp/%s/README' % buildDir, tDownload,
                   cdir='software/tmp')
    SideEffect('dummy', tUntar)  # so it works fine in parallel builds
    Clean(tUntar, 'software/tmp/%s' % buildDir)
    lastTarget = tUntar

    for target, command in extraActions:
        lastTarget = env.Command(
            target,
            lastTarget,
            Action(command, chdir='software/tmp/%s' % buildDir))
        SideEffect('dummy', lastTarget)  # so it works fine in parallel builds

    # Clean the special generated files.
    for cFile in clean:
        Clean(lastTarget, cFile)

    # Add the dependencies.
    for dep in deps:
        Depends(tUntar, dep)

    if default or GetOption(name):
        Default(lastTarget)

    return lastTarget


# Add methods so SConscript can call them.
env.AddMethod(addLibrary, "AddLibrary")
env.AddMethod(addModule, "AddModule")
env.AddMethod(addPackage, "AddPackage")
env.AddMethod(manualInstall, "ManualInstall")


#  ************************************************************************
#  *                                                                      *
#  *                            Extra options                             *
#  *                                                                      *
#  ************************************************************************


opts = Variables(None, ARGUMENTS)

opts.Add('MPI_CC', 'MPI C compiler', 'mpicc')
opts.Add('MPI_CXX', 'MPI C++ compiler', 'mpiCC')
opts.Add('MPI_LINKERFORPROGRAMS', 'MPI Linker for programs', 'mpiCC')
opts.Add('MPI_INCLUDE', 'MPI headers dir ', '/usr/include')
opts.Add('MPI_LIBDIR', 'MPI libraries dir ', '/usr/lib')
opts.Add('MPI_LIB', 'MPI library', 'mpi')
opts.Add('MPI_BINDIR', 'MPI binaries', '/usr/bin')

opts.Add('SCIPION_HOME', 'Scipion base directory', abspath('.'))

opts.Update(env)

# TODO: check the code below to see if we can do a nice "purge".

# def _removeInstallation(env):
#     """
#     Function that cleans the folders used by a scipion installation in order to completely remove everything related to that installation
#     """
#     # Dictionary to store the folder that need to be emptied (TOCLEAN) or deleted (TOREMOVE)
#     UNINSTALL = {'TOCLEAN': [join('software','lib'),
#                              join('software', 'lib64'),
#                              join('software', 'bin'),
#                              join('software', 'man'),
#                              join('software', 'share'),
#                              join('software', 'tmp'),
#                              join('software', 'log')],
#                  'TOREMOVE': [join('software', 'install', 'scons-2.3.1')]}
# #    if _ask("Proceeding with Scipion purge process. Everything is going to be removed from the machine. Are you sure?") != 'y':
# #        return False
#     for dir in UNINSTALL.get('TOCLEAN'):
#         print "Cleaning %s" % dir
#         list = os.listdir(dir)
#         for thing in list:
#             path = join(dir, thing)
#             if thing == '.gitignore':
#                 continue
#             if os.path.isfile(path) or os.path.islink(path):
#                 os.unlink(path)
#             else:
#                 shutil.rmtree(path)
#     for dir in UNINSTALL.get('TOREMOVE'):
#         print "Deleting %s" % dir
#         shutil.rmtree(dir)
#     return True

# ########################
# # Command-line options #
# ########################

# AddOption('--update',
#           dest='update',
#           action='store_true',
#           help='Check for packages or libraries updates')
# AddOption('--purge',
#           dest='purge',
#           action='store_true',
#           help='Completely clean the installation and its binaries')
# AddOption('--binary',
#           dest='binary',
#           action='store_true',
#           help='After doing the installation, create and package a binary for distribution')

AddOption('--with-all-packages', dest='withAllPackages', action='store_true',
          help='Get all EM packages')


Export('env')

env.SConscript('SConscript')
