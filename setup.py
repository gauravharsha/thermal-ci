from numpy.distutils.core import setup, Extension


thermalcisd = Extension(
    'tfdcisd.ThermalCISD',
    sources=['tfdcisd/fort_src/ThermalCISD.f90', ],
    extra_f90_compile_args=["-O3", "-lblas", "-llapack"]
    # extra_f90_compile_args=["-O3", "-fopenmp", "-lgomp", "-lblas"]
)

expvals = Extension(
    'tfdcisd.ExpVals',
    sources=['tfdcisd/fort_src/ExpVals.f90', ],
    extra_f90_compile_args=["-O3", "-lblas", "-llapack"]
    # extra_f90_compile_args=["-O3", "-fopenmp", "-lgomp", "-lblas"]
)

setup(
   name='tfdcisd',
   version='0.1.1',
   description='A package for calculating finite temperature configuration\
       interaction theory in the grand canonical ensemble',
   packages=['tfdcisd'],
   ext_modules=[thermalcisd, expvals]
)
