# How to install MSS for MATLAB

The Marine Systems Simulator (MSS) is a MATLAB Toolbox for for designing and testing marine control systems. The M-files are also compatible with [GNU Octave](https://www.octave.org). 

Choose one of the three installation methods below:

## 1. Install from inside MATLAB (recommended)

1. Open the **Apps** tab in MATLAB.
2. Click **Get More Apps**.
3. Search for `MSS`.
4. Select **Marine Systems Simulator (MSS)** and click **Install**.

## 2. Download the ZIP archive from GitHub

1. [Download the MSS ZIP archive](https://github.com/cybergalactic/MSS/archive/refs/heads/master.zip).
2. Extract the archive to a desired location and name the directory `MSS`.

## 3. Install from the MATLAB Central webpage

1. Open [Marine Systems Simulator (MSS) on MATLAB Central](https://www.mathworks.com/matlabcentral/fileexchange/86393-marine-systems-simulator-mss).
2. Click **Download** or **Install in MATLAB** and follow the prompts.

For MATLAB R2026b or later, the package can also be installed from the MATLAB Command Window:

```matlab
mpminstall("marine_systems")
```

## Set up the MATLAB path

After installing or downloading MSS, run:

```matlab
mssPath
```

If MATLAB cannot find `mssPath`, add MSS to the path first:

1. On the MATLAB **Home** tab, select **Set Path**.
2. Select **Add with Subfolders**, choose the MSS directory, and save the path.

Alternatively, open the MSS directory as the current folder and run:

```matlab
addpath(genpath(pwd))
savepath
mssPath
```

The `mssPath` command refreshes the MSS folders on the MATLAB path, saves the path, and removes obsolete MSS path entries. Run it again after updating MSS.

## Get started

Display the MSS help menu from the MATLAB Command Window:

```matlab
mssHelp
```
