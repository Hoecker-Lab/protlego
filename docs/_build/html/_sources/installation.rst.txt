Installation
============

Protlego can be installed via github, pip, or conda.

PIP nstallation
===============

To install Protlego via pip, type:
 .. code-block:: javascript

    pip install protlego

Graph-tool is not included in Pypi, you will need to install it separately. To do so, please refer to the graph-tool documentation: https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions



CONDA installation
==================

To install Protlego via conda, you will first have to install conda. If you are familiar with conda you might already have an environment, otherwise create one as follows:

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

You will need to adjust the boost installation after installing. Type:

 .. code-block:: javascript

	conda install -c anaconda boost=1.69

Also make sure your cairo package comes from channel conda-forge. If it is not, you can type:

 .. code-block:: javascript

        conda install -c conda-forge cairo

Once the environment installed, you can start designing chimeras with protlego!

Setup
============

You will need a Fuzzle database. Download the latest version from here:
https://132.180.65.101/static/js/hh207clusters.csv

Then go to your protlego folder, and run the script:
 .. code-block:: javascript

       python database/database_setup.py

Make sure you have changed the path to your hh207clusters.csv accordingly.
