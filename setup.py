#!/usr/bin/env python
'''
MAGeCK2 set up script
'''


from __future__ import print_function;

import os
import sys
from distutils.core import setup, Extension
from subprocess import call as subpcall
from distutils.command.install import install as DistutilsInstall

try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py

exec(open('mageck2/version.py').read())

def compile_rra():
  #
  os.chdir('rra');
  subpcall('make',shell=True);
  rev=subpcall('bin/RRA',shell=True);
  os.chdir('../');
  return rev;


class RRAInstall(DistutilsInstall):
  def run(self):
    # compile RRA
    if(compile_rra()!=0):
      print("CRITICAL: error compiling the RRA source code. Please check your c compilation environment.",file=sys.stderr);
      sys.exit(1);
    DistutilsInstall.run(self)



def main():
  # check python version
  if float(sys.version[:3])<2.7:
    sys.stderr.write("CRITICAL: Python version must be >=2.7!\n")
    sys.exit(1);

  setup(name='mageck2',
    version=__version__,
    description='Model-based Analysis of Genome-wide CRISPR-Cas9 Knockout',
    author='Wei Li ',
    author_email='li.david.wei@gmail.com',
    url='https://github.com/davidliwei/mageck2',
    packages=['mageck2'],
    scripts=['bin/mageck2'],
    package_dir={'mageck2':'mageck2'},
    cmdclass={'install':RRAInstall, 'build_py': build_py},
    package_data={'mageck2':['*.Rnw','*.RTemplate','*.gmt','*.txt','*.Rmd']},
    data_files=[('bin', ['rra/bin/RRA'])]
    #package_data={'mageck':['mageck/Makefile','mageck/src/*.c','include/*','utils/*']}
    #data_files=[('',['Makefile','src/*.c','include/*','utils/*'])]
  );


if __name__ == '__main__':
  main();
