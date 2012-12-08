from distutils.core import setup, Extension
import numpy as np

setup(
    name='Nurbs',
    version='0.1',
    description='Python module to work with NURBS curves and surfaces.',
    author='Runar Tenfjord',
    author_email='runten@netcom.no',
    url='http://runten.tripod.com/',
    packages=['Nurbs', 'Nurbs.demos'],
    include_dirs = [np.get_include()],
    ext_modules = [Extension("Nurbs._Bas", ["Nurbs/_Bas.c"])],
    data_files=[('Nurbs/Doc', ['LICENSE', 'README'])]
    )

