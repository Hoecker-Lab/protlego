Installation
============

Protlego can be installed via github, pip, or conda.

PIP installation
================

To install Protlego via pip, type:
 .. code-block:: javascript

    pip install protlego

Graph-tool is not included in Pypi, you will need to install it separately. To do so, please refer to the graph-tool documentation: https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions



CONDA installation
==================

We only support Linux at the moment. To install Protlego via conda, you will first have to install conda. If you are familiar with conda you might already have an environment, otherwise create one as follows:

 .. code-block:: javascript

    conda create -n myenv

And then activate it:

 .. code-block:: javascript

    conda activate myenv

Later, you might need to add some channels:

 .. code-block:: javascript

    conda config --append channels bioconda
    conda config --append channels conda-forge
    conda config --append channels acellera

And then you are good to install protlego:

 .. code-block:: javascript

    conda install -c nferruz protlego


Once the environment installed, you can start designing chimeras with protlego!

Setup
============

You will need a Fuzzle database. You can download the latest version here: 
https://fuzzle.uni-bayreuth.de:8443/static/fuzzle2.07.db


Then place it in your protle place it in your protlego folder. If you don't know where
your protlego folder got installed, import protlego
in a python terminal and see its path:

 .. code-block:: bash

 	$ python
	>>> import protlego
	>>> protlego.__file__


Installing from environment file
=================================
If you have issues installing Protlego, please post an issue in our GitHub repository and we will assist you.
Alternatively, you can create a conda environment directly from the requirements. 

Download the environment file from our repository: https://github.com/Hoecker-Lab/protlego
And then:

 .. code-block:: javascript

    conda env create -f environment.yml


Troubleshooting
===============
Post issues in our GitHub repository:
https://github.com/Hoecker-Lab/protlego/issues/


