# How to install MSS for MATLAB

The Marine Systems Simulator (MSS) is a Matlab and Simulink library for marine systems (https://www.mathworks.com).

 Installation
-
Download the MSS Directory:
1. Navigate to the GitHub MSS repository and click the green "Code" button.
2. Select "Download ZIP" from the dropdown menu.
3. Unzip the file to a desired location on your computer once the download is complete.

Set Up MATLAB Path:
1. Open MATLAB.
2. From the MATLAB menu, choose "Set Path."
3. Select "Add with Subfolders" and navigate to the unzipped MSS directory. Alternatively, suppose you are already in the top-level folder of the MSS toolbox. In that case, you can set the path directly by executing the following command in the MATLAB command window:

       addpath(genpath(pwd))
       savepath

Later, you can update an existing path automatically and remove dead links by using the command

    mssPath