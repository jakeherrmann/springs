
DESCRIPTION:

This repository contains code for (1) solving equilibrium configurations of elastic spring networks, and (2) simulating agent-based models of respiratory biomechanics.  The core spring network solver is written in C++ and must be compiled.  There are two interfaces: python and matlab.  Both contain functionality to compile and execute the C++ solver, with data communication via temporary input and output files.  Data post-processing and visualization is also available through the matlab and python interfaces.

This is very much a work in progress.  Any requests to use this code for other purposes should be directed to Jacob Herrmann ( jakeherr @ bu . edu ) and Bela Suki ( bsuki @ bu . edu ).


SETUP:

The Python and Matlab commands to setup the solver executable (springs_setup.py and springs_setup.m, respectively) expect the GNU make and g++ commands to be accessible from the command prompt on your operating system.

	- Windows OS: We recommend MSYS and MinGW (https://www.msys2.org/).  Be sure to add the msys64\usr\bin and msys64\mingw64\bin folders to your PATH environment variable.

	- Mac OS: Users may install XCode or just the command line developer tools (http://www.edparrish.net/common/macgpp.php).


USING GIT:

Pulling new updates from the git repository will be much easier if you do not make any changes to the files provided.  Any errors or bugs in the master files should be corrected by someone with write privileges to the master respository.  The git ignore is set to ignore any files or folders that match the name “test_*” (for testing and debugging), “dev_*” (for development and documentation), or “user_*” (for personal projects).  I recommend you to use these prefixes on names of files and folders that are intended for your personal use.