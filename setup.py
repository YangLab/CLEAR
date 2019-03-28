#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from src.version import __version__

setup(name='CLEAR',
      version=__version__,
      description='Circular RNA quantification toolkits',
      author='Xu-Kai Ma',
      author_email='maxukai@picb.ac.cn',
      url='https://github.com/YangLab/CLEAR',
      license='MIT',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
      ],
      keywords='circular RNAs',
      packages=find_packages(),
      install_requires=[
          'pysam>=0.8.4',
          'pybedtools',
      ],
      entry_points={
          'console_scripts': [
              'fast_quant=scr.run:main',
              'circ_quant=scr.circ_quant:main',
          ],
      },
      )
