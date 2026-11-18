# Configure WAMIT for MSS import

The bundled tanker is a worked example for a single vessel at zero
forward speed with port–starboard symmetry.

## ASCII input files

Copy [`fnames.wam`](tanker_wamitData/fnames.wam),
[`tanker.cfg`](tanker_wamitData/tanker.cfg),
[`tanker.pot`](tanker_wamitData/tanker.pot),
[`tanker.frc`](tanker_wamitData/tanker.frc), and
[`tanker.gdf`](tanker_wamitData/tanker.gdf) to one working folder on the PC.
On GitHub, open a file and use **Download raw file** to save it with the
exact name. The blocks below show the WAMIT 6 inputs used for the bundled tanker results.

### [`fnames.wam`](tanker_wamitData/fnames.wam)

```text
tanker.cfg
tanker.pot
tanker.frc
tanker.gdf
```

WAMIT reads this list to find the case files. Keep the listed filenames
with their original extensions.

### [`tanker.cfg`](tanker_wamitData/tanker.cfg)

```text
IALTFRC=2
IALTPOT=1
ISOLVE=1
MAXSCR=5000
ILOWHI=0
MAXITT=100
IPERIO=1
USERID_PATH=c:\WAMITv63
```

Set `USERID_PATH` to the directory containing your WAMIT license file.
`IALTPOT=1` and `IALTFRC=2` identify the layouts shown below. WAMIT 6 uses
`IPERIO=1` for periods in seconds. For WAMIT 7, convert these WAMIT 6
inputs with the supplied `v6v7inp` utility and check the converted files.
Use `IPERIN=1` and `IPEROUT=1` to preserve periods in seconds and
`ISOR=1` for low-order pressure integration (`.9`).
See Chapters 4 and 6 and Appendix B of the
[WAMIT User Manual](https://www.wamit.com/manual7.x/v75_manual.pdf).

### [`tanker.pot`](tanker_wamitData/tanker.pot)

```text
Tanker Lpp = 246 m, B = 46 m, T = 10 m
1
3
-1 0.00 0.00 0.00 0.00
0 1 0
1 1
1 1 1 1 1 1
60
-1 0 1.05 1.10 1.17 1.24 1.32 1.40 1.50 1.60 1.72 1.85 1.99 2.14 2.32 2.51 2.72 2.95 3.21 3.49 3.81 4.16 4.54 4.96 5.43 5.95 6.53 7.16 7.86 8.63 9.49 10.44 11.48 12.64 13.91 15.32 16.88 18.60 20.51 22.61 24.94 27.51 30.35 33.49 36.95 40.79 45.03 49.71 54.88 60.60 66.92 73.91 81.63 90.16 99.59 110.01 121.53 134.26 148.33 163.88
19
0   10   20   30   40   50   60   70   80   90  100  110  120  130  140  150  160  170  180
```

The six ones enable the rigid-body DOFs. The list of 60 periods begins
with `-1` (infinite period) and `0` (zero period). The remaining periods
are in seconds. The downloadable file also lists two trailing periods
(`181.06` and `200.05` s). The block ends at `163.88` s to match the 60
entries used in the bundled outputs. Keep the heading count `19` on a
line by itself with no
trailing spaces, followed by the 19 headings from 0° to 180° in 10° steps.
MSS uses port–starboard symmetry to fill in the other headings.

### [`tanker.frc`](tanker_wamitData/tanker.frc)

```text
Tanker Lpp = 246 m, B = 46 m, T = 10 m
1 0 1 1 0 0 0 0 1
1025.00
3.93 0.00 2.5
1
94620210.00 0.00 0.00 0.00 236550525.00   0.00
0.00 94620210.00 0.00 -236550525.00 0.00   0.00
0.00 0.00 94620210.00 0.00 -371857425.30   0.00
0.00 -236550525.00 0.00 27409620280.88 0.00   0.00
236550525.00 0.00 -371857425.30 0.00 387080076077.14   0.00
0.00 0.00 0.00 0.00 0.00 417428070207.44
1
410000.00 0.00 0.00 0.00 0.00 0.00
0.00 0.00 0.00 0.00 0.00 0.00
0.00 0.00 0.00 0.00 0.00 0.00
0.00 0.00 0.00 1560000000.00 0.00 0.00
0.00 0.00 0.00 0.00 0.00 0.00
0.00 0.00 0.00 0.00 0.00 0.00
1
0.00 0.00 0.00 0.00 0.00 0.00
0.00 0.00 0.00 0.00 0.00 0.00
0.00 0.00 0.00 0.00 0.00 0.00
0.00 0.00 0.00 0.00 0.00 0.00
0.00 0.00 0.00 0.00 0.00 0.00
0.00 0.00 0.00 0.00 0.00 0.00
0
0
```

The option line requests added mass and damping (`.1`), diffraction
force RAOs (`.3`), motion RAOs (`.4`), and diagonal pressure mean drift
(`.9`). The next two lines give water density in kg/m³ and the center of
gravity in metres. The first 6-by-6 matrix is the rigid-body mass matrix
about the WAMIT origin. The FPSO uses
[`fpso.frc`](fpso_wamitData/fpso.frc) with alternative 1, which provides
the vertical center of gravity and radii of gyration.

### [`tanker.gdf`](tanker_wamitData/tanker.gdf)

```text
Tanker Lpp = 246 m, B = 46 m, T = 10 m
1 9.80665 	ULEN GRAV
0  1 	ISX  ISY
2375
```

Only the geometry header is shown here. Download the complete file for
the 2,375 panels. The WAMIT axes are forward, port, and up, with the
origin at midships on the centerline in the still-water plane.
MSS converts the axes to forward, starboard, and down and retains the
WAMIT origin for the wave loads.

## Run WAMIT

Open a command prompt in the working folder and run `wamit` using the
installed executable. The bundled case uses `fnames.wam` to select the
inputs. Keep the original MSS copy of the case as a reference when running
WAMIT on a PC.

## Files to keep for import

After running WAMIT, keep the following ASCII files from the same run together.
The bundled case uses the base name `tanker`.

| Files | Data read by `wamit2vessel` |
| --- | --- |
| `.pot`, `.frc`, `.out` | Heading count, mass properties, and formatted results including hydrostatics and length scale |
| `.1` | Added mass and radiation damping, including the `-1` and `0` period entries |
| `.3` or `.2` | Excitation force RAOs from diffraction or Haskind |
| `.4` | Motion RAOs |
| `.9` or `.8` | Mean drift from pressure or momentum integration. Use diagonal output (`IOPTN(9)=1` or `IOPTN(8)=1`) with surge, sway, and yaw at every positive period and heading in the drift grid. |

## Import into MSS

Add MSS to the MATLAB or GNU Octave path. Set `caseFolder` to the folder
containing the run files. The call takes the shared filename without an
extension, followed by draught, length between perpendiculars, and overall
breadth in metres.

```matlab
caseFolder = 'C:\WAMIT\cases\tanker';
vessel = wamit2vessel(fullfile(caseFolder,'tanker'),10,246,46,'0000', ...
    [49.6 3.6 4.5],[0 0.1 0],[5 0 0]);
```

Press Return at the alternative 2 mass check. The importer saves
`tanker.mat` beside `wamit2vessel.m`, replacing any existing file with
that name. If both supported force or drift files are present, the importer
asks which one to use. For a new vessel, supply or enter vessel-specific
surge, sway, and yaw decay time constants. With FRC alternative 1, you can
supply water density in kg/m³ as the ninth argument.

See [WAMIT vessel data](WAMIT-Vessels.md) for the other bundled import calls
and the damping arguments.
