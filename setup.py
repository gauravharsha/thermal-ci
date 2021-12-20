from numpy.distutils.core import setup, Extension


thermalcisd = Extension(
    'tfdcisd.ThermalCISD',
    sources=['tfdcisd/fort_src/ThermalCISD.f90', ],
    libraries=['lapack', 'blas'],
    extra_f90_compile_args=["-O3"]
    # extra_f90_compile_args=["-O3", "-fopenmp", "-lgomp", "-lblas"]
)

expvals = Extension(
    'tfdcisd.ExpVals',
    sources=['tfdcisd/fort_src/ExpVals.f90', ],
    libraries=['lapack', 'blas'],
    extra_f90_compile_args=["-O3"]
    # extra_f90_compile_args=["-O3", "-fopenmp", "-lgomp", "-lblas"]
)

fixrefthermalcisd = Extension(
    'tfdcisd.FixRefThermalCISD',
    sources=['tfdcisd/fort_src/FixRefThermalCISD.f90', ],
    extra_f90_compile_args=["-O3"]
)

fixrefexpvals = Extension(
    'tfdcisd.FixRefExpVals',
    sources=['tfdcisd/fort_src/FixRefExpVals.f90', ],
    extra_f90_compile_args=["-O3"]
)

setup(
   name='tfdcisd',
   version='1.2',
   description='A package for calculating finite temperature configuration\
       interaction theory in the grand canonical ensemble',
   packages=['tfdcisd'],
   ext_modules=[thermalcisd, expvals, fixrefthermalcisd, fixrefexpvals]
)
