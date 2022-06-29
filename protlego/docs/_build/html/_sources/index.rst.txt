.. ProtLego documentation master file, created by
   sphinx-quickstart on Thu Jun  6 19:37:27 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ProtLego! An Open-source Python Library for Chimera Design and Analysis
==================================================================================
ProtLego is a tool that enables constructing protein chimeras and its structural analysis.
ProtLego is writen in Python, so it facilitates that other scientists can extend it to their needs.

With ProtLego, it is possible to go from a pair of PDBs that share a structural fragment in common, to build
all possible chimeras between the two, rank them with the AMBER or CHARMM molecular mechanics forcefields,
and analyze their structural features (hydrogen networks, salt bridges, hydrophobic clusters...etc.)

.. toctree::
   :maxdepth: 2
   :caption: Code 
   
   builder
   database
   networks
   energy
   analyse


.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   fuzzle

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   introduction
   fetching
   network
   building
   energy
   structural


Publications
============

**Protlego: A Python package for the analysis and design of chimeric proteins.** 
Noelia Ferruz, Jakob Noske, Birte Höcker.
