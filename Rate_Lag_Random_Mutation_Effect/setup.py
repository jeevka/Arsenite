from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("Yeast_Simulator", ["Yeast_Simulator.pyx"], include_dirs =[ np. get_include ()])]
)

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("Yeast_Simulator_Cell_Division", ["Yeast_Simulator_Cell_Division.pyx"], include_dirs =[ np. get_include ()])]
)


setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("Yeast_Simulator_Subprograms", ["Yeast_Simulator_Subprograms.pyx"], include_dirs =[ np. get_include ()])]
)

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("Yeast_Simulator_Mutations", ["Yeast_Simulator_Mutations.pyx"], include_dirs =[ np. get_include ()])]
)


#setup(
#    cmdclass = {'build_ext': build_ext},
#    ext_modules = [Extension("Yeast_Simulator_Gene_Duplications", ["Yeast_Simulator_Gene_Duplications.pyx"], include_dirs =[ np. get_include ()])]
#)
