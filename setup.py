from distutils.core import setup
from setuptools import find_packages

setup(
  name = 'protlego',
  packages = find_packages(),
  include_package_data=True,
  version = '1.0.12',
  description = 'A python package for the automated construction of chimeras and analysis',
  author = 'Noelia Ferruz',                   # Type in your name
  author_email = 'noelia.feruz@gmail.com',      # Type in your E-Mail
  url = 'https://github.com/Hoecker-Lab/protlego',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/Hoecker-Lab/protlego/archive/1.0.10.tar.gz',
  keywords = ['PROTEIN DESIGN', 'CHIMERAGENESIS', 'PROTEIN EVOLUTION', 'NETWORK THEORY', 'FUZZLE'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
            'numpy',
            'scipy',
            'csb',
            'moleculekit',
            'matplotlib',
            'openmm',
            'mdtraj',
            'pandas',
            'boost',
    ],
  classifiers=[
    'Development Status :: 5 - Production/Stable',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
	'Programming Language :: Python :: 3.6',
	'Programming Language :: Python :: 3.7',
	'Programming Language :: Python :: 3.8',
	'Programming Language :: Python :: 3.9',
  ],
)
