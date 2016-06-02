from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        Extension("integrals",
                  ["integrals.pyx"],
                  ),
    ])
setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        Extension("files_operations",
                  ["files_operations.pyx"],
                  ),
    ])
setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        Extension("integrals_analytical",
                  ["integrals_analytical.pyx"],
                  ),
    ])
setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        Extension("orthogonalization",
                  ["ortogonalization.pyx"],
                  ),
    ])