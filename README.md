# VirusDecode Bioinformatics

This repository implements bioinformatics functionalities for the VirusDecode project, focusing on sequence alignment using Python and MUSCLE. We also use opensource LinearDesign as well.

## Development Environment

### Operating System
- Linux

### Programming Language
- Python 3.11.9
- Python 2.7 for running LinearDesign

### Required Packages
- Biopython 1.83

### Alignment Tool
- MUSCLE v3.8.1551 by Robert C. Edgar

## Installation Instructions
To install the required biopython and muscle, run the following command:
```sh
sudo apt update
sudo apt install python2
sudo apt install python3
pip install biopython

mkdir muscle
cd muscle
wget https://www.drive5.com/muscle/muscle_src_3.8.1551.tar.gz
tar -xvzf muscle_src_3.8.1551.tar.gz
make
sudo cp muscle /usr/local/bin/
cd ..
rm -rf muscle
```

## Execution
```sh
git clone https://github.com/LinearDesignSoftware/LinearDesign.git
cd LinearDesign
make
cd..
git clone https://github.com/VirusDecode-project/bioinformatics.git
cd bioinformatics
python3 virusdecode.py
```
