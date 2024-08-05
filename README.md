# VirusDecode Bioinformatics

This repository implements bioinformatics functionalities for the VirusDecode project, focusing on sequence alignment using Python and MUSCLE.

## Development Environment

### Operating System
- Linux

### Programming Language
- Python 3.11.9

### Required Packages
- Biopython 1.83

### Alignment Tool
- MUSCLE v3.8.1551 by Robert C. Edgar

## Installation Instructions
To install the required biopython and muscle, run the following command:
```sh
pip install biopython
sudo apt-get update

mkdir muscle
cd muscle
wget https://www.drive5.com/muscle/muscle_src_3.8.1551.tar.gz
tar -xvzf muscle_src_3.8.1551.tar.gz
make
sudo cp muscle /usr/local/bin/
```
