# protlego

Installation
============

Protlego can be installed via github, pip, or conda.

PIP installation
----------------

To install protlego via pip, just type the following command:

`pip install protlego`

Taking into account that graph-tool is not included in Pypi, you will need to install it separately. To do so, please refer to the graph-tool documentation: https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions

CONDA installation
----------------

To install protlego via conda, type the following:

`conda install -c anaconda protlego`

You will need to adjust the boost installation after installing. Type:

`conda install -c anaconda boost=1.69`

Also make sure your cairo package comes from channel conda-forge. If it is not, you can type:

`conda install -c conda-forge cairo`

Once the environment installed, you can start designing chimeras with protlego!

Setup
============

You will need a Fuzzle database. Download the latest version from here: 

https://132.180.65.101/static/js/hh207clusters.csv.

Then go to your protlego folder, and run the script 

`python database/database_setup.py`. 

Make sure you have changed the path to your hh207clusters.csv accordingly.
