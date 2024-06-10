# How to install MSS for GNU Octave

The m-files of the MSS (Marine Systems Simulator) toolbox are compatible with MATLAB and the free software GNU Octave (https://www.octave.org), facilitating broad accessibility and application in marine systems simulation. To begin using the MSS toolbox, ensure it is installed and properly set up in your GNU Octave environment. 

## Download the MSS Directory:
1. [Download ZIP](https://github.com/cybergalactic/MSS/archive/refs/heads/master.zip)
2. Unzip the file to a desired location on your computer and name the directory MSS.
   
## Download and Install GNU Octave:
1. Download and install GNU Octave from https://octave.org/download.
2. Edit the start-up file (https://docs.octave.org/latest/Startup-Files.html) to include

       graphics_toolkit("qt");

   The function call 'graphics_toolkit ("qt")' selects the Qt/OpenGL system.
3. Start GNU Octave from the terminal using 'octave --gui'. Type

       >> graphics_toolkit
       ans = qt

   to verify that you are using QT.
4. Select "Edit - Set path - Add directory with subdirectories" and navigate to the unzipped MSS directory. Select "Open - Save - Close".
5. To get started and find help on using the MSS, type the following command in the command window:
   
       >> mssHelp
