#! /usr/bin

# @Author: Yu Zhu
# @Email: 1830416012@stu.suda.edu.cn
# @Address: Department of Bioinformatics, Medical College, Soochow University

## Special Attention: All the installation procedures may cause various errors,
##					  Please check feedback and make adjustments to the following installation commands.

## ===================================================================================================== ##

# Prepare installation tools ready
sudo apt-get update
sudo apt-get upgrade

# Install required files/tools/softwares
sudo apt-get install cmake

## @ Attention: 
##			apt-get command may not install cmake with proper version(>= 3.20.0), 
##			alternative installation method is required as follows：
##
##			  @ Remove installed cmake version:
## 				sudo apt-get autoremove cmake
## 				cd
##
##			  @ Download cmake packages in .tar.gz format: 
## 				wget https://cmake.org/files/v3.20/cmake-3.20.0.tar.gz
## 				tar xvf cmake-3.20.0.tar.gz
## 				cd cmake-3.20.0
##
##			  @ Make installation:
## 				./bootstrap --prefix=/usr
## 				make
## 				sudo make install
##
##			  @ Check the current cmake version: (For latest Gromacs, the cmake version must exceed 3.20.0)
## 				cmake --version

# Install required files/tools/softwares
sudo apt-get install build-essential

# Make new directory for downloading
mkdir Gromacs
cd Gromacs

# Download Gromacs packages and regression test files in .tar.gz format:
wget ftp://ftp.gromacs.org/gromacs/gromacs-2021.1.tar.gz
wget https://ftp.gromacs.org/regressiontests/regressiontests-2021.1.tar.gz
tar xvzf gromacs-2021.1.tar.gz
tar xvzf regressiontests-2021.1.tar.gz

# Install required files/tools/softwares
sudo apt-get install libfftw3-dev

# Make new directory for installation
cd gromacs-2021.1/
mkdir build
cd build

# Install Gromacs
sudo cmake .. -DGMX_BUILD_OWN_FFTW=OFF -DREGRESSIONTEST_DOWNLOAD=OFF -DCMAKE_C_COMPILER=gcc -DREGRESSIONTEST_PATH=/root/Gromacs/regressiontests-2021.1
sudo make check
sudo make install
source /usr/local/gromacs/bin/GMXRC

# Check if installation has been done
gmx
gmx pdb2gmx --version

## @ Attention: 
##			In order to execute source command every time before using Gromacs,
##			modifying environment PATH is recommended.
##			Related commands are shown below:
##		
##			  @ Change working directory into a default one:
## 				cd
##
##			  @ Modify PATH in .bashrc file with vim tool :
## 				vim .bashrc
##
##			  @ Append the export line in .bashrc file and quit:
## 				export PATH=$PATH:/usr/local/gromacs/bin
## 				
##			  @ Refresh the current shell environment:
##				source .bashrc

# Check if installation has been done
gmx
gmx pdb2gmx --version

