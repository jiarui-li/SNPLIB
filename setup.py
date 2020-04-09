import os
import re
import sys
import platform
import subprocess

from setuptools import setup, find_packages, Command, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


def get_extra_cmake_options():
    _cmake_extra_options = []

    opt_key = None
    argv = [arg for arg in sys.argv]

    for arg in argv:
        if opt_key == 'use_mkl':
            _cmake_extra_options.append('-DUSE_MKL')
        elif opt_key == 'intel_root':
            _cmake_extra_options.append(
                '-DINTEL_ROOT={arg}'.format(arg=arg.strip()))
        elif opt_key == 'use_openblas':
            _cmake_extra_options.append('-DUSE_OPENBLAS')
        elif opt_key == 'openblas_root':
            _cmake_extra_options.append(
                '-DOPENBLAS_ROOT={arg}'.format(arg=arg.strip()))

        if opt_key:
            sys.argv.remove(arg)
            opt_key = None
            continue

        if arg in ['--use_mkl', '--use_openblas', '--intel_root', '--openblas_root']:
            opt_key = arg[2:].lower()
            sys.argv.remove(arg)
            continue

    return _cmake_extra_options


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(
            self.get_ext_fullpath(ext.name)))
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cmake_args += get_extra_cmake_options()

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += [
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] +
                              cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] +
                              build_args, cwd=self.build_temp)


setup(
    name="SNPLIB",
    license="BSD 3-Clause License",
    packages=['SNPLIB'],
    ext_modules=[CMakeExtension('_SNPLIB')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
