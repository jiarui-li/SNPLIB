import os
import re
import sys
import platform
import subprocess

from setuptools import setup, find_packages, Command, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


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

        blas_library = input('Use Intel MKL(m) or OpenBLAS(o)?')
        blas_library = blas_library.lower()[0]
        while blas_library not in ['m', 'o']:
            blas_library = input('Use Intel MKL(m) or OpenBLAS(o)?')
            blas_library = blas_library.lower()[0]

        pf = platform.system()
        intel_root = None
        if blas_library == 'm':
            cmake_args += ['-DUSE_MKL=ON']
            if pf == "Windows":
                intel_default_root = 'C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows'
            if pf == "Linux":
                intel_default_root = '/opt/intel'
            if pf == "Darwin":
                intel_default_root = '/opt/intel'
            intel_root = input(
                'The directory of Intel MKL [{}]: '.format(intel_default_root))
            if not intel_root:
                intel_root = intel_default_root

            if not os.path.exists(intel_root):
                raise RuntimeError(
                    "Intel MKL must be installed in {}".format(intel_root))
            cmake_args += ['-DINTEL_ROOT={}'.format(intel_root)]
        elif blas_library == 'o':
            cmake_args += ['-DUSE_OPENBLAS=ON']
            if pf == "Windows":
                openblas_default_root = 'C:/opt'
            if pf == "Linux":
                openblas_default_root = '/opt/OpenBLAS'
            if pf == 'Darwin':
                openblas_default_root = '/opt/OpenBLAS'
            openblas_root = input(
                'The directory of OpenBLAS [{}]: '.format(openblas_default_root))
            if not openblas_root:
                openblas_root = openblas_default_root

            if not os.path.exists(openblas_root):
                raise RuntimeError(
                    "OpenBLAS must be installed in {}".format(openblas_root))
            cmake_args += ['-DOPENBLAS_ROOT={}'.format(openblas_root)]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]
        if pf == "Windows":
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
    version="0.1.0",
    author="Li Jiarui",
    author_email="jiarui_li@outlook.com",
    packages=['SNPLIB'],
    license="BSD 3-Clause License",
    ext_modules=[CMakeExtension('_SNPLIB')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    install_requires=["numpy", "scipy", "pandas"],
)
