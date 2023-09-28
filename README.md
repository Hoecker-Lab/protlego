# protlego

Protlego is a evolution based chiemragenesis tool. 
It works by identifying shared remote homology between two unrelated domains and generates chimeras by mimicking eith of the two evolutionary event :
Recombination and Insertion 

After generatin gthe chimeras, user can downstream minimise the models and also perform various structural analysis.

Installation
============

Follow these steps to install and run the tool on your local machine.

1. Clone the github repository :

  ```bash
  git clone -n add-features https://github.com/Hoecker-Lab/protlego/tree/add_features
  ```

2. Download fuzzle database using the following link below :

  https://fuzzle.uni-bayreuth.de:8443/static/fuzzle2.07.db

3. cd into the folder where protlego is installed and move the fuzzle database file to ./protlego/database folder.

4. Create conda environment using the command :

 ```bash
 conda env create -f simple_env.yml
 ```

5. Activate conda environment 

 ```bash
 conda activate protlego
 ```

6. Follow the latest guidelines at: https://hoecker-lab.github.io/protlego/installation.html
