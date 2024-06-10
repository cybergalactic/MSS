# How to install MSS for MATLAB

The Marine Systems Simulator (MSS) is a Matlab and Simulink library for marine systems (https://www.mathworks.com). The m-files are compatible with the free software GNU Octave (https://www.octave.org). To begin using the MSS toolbox, ensure it is installed and properly set up in your MATLAB environment.

## Download the MSS Directory:
1. [Download ZIP](https://github.com/cybergalactic/MSS/archive/refs/heads/master.zip)
2. Unzip the file to a desired location on your computer and name the directory MSS.

## Set Up MATLAB Path:
1. Open MATLAB.
2. From the MATLAB menu, choose "Set Path."
3. Select "Add with Subfolders" and navigate to the unzipped MSS directory. Alternatively, suppose you are already in the top-level folder of the MSS toolbox. In that case, you can set the path directly by executing the following command in the MATLAB command window:

       >> addpath(genpath(pwd))
       >> savepath

Later, you can update an existing path automatically and remove dead links by using the command

       >> mssPath
    
4. To get started and find help on using the MSS, type the following command in the command window:

       >> mssHelp
