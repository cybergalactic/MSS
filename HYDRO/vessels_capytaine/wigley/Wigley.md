# Wigley benchmark hull

`wigley.mat` is copied from `MSS-Capytaine/vessels_capytaine/wigley_stl/wigley.mat`.
It contains zero-speed Capytaine hydrodynamics, first-order RAOs, signed mean
drift and JSON-selected damping inputs under `vessel.powerBased`.
The saved inputs are `T_126 = [80 8 12]` s,
`delta_zeta_345 = [0 0.1 0]` and `nonlinear_456 = [0 0 0]`.

Select the file through **Load vessel MAT-file** in `SIMhydroVessel`.
Regenerate with `python main.py wigley_stl` in MSS-Capytaine, or refresh only
damping inputs with `python main.py wigley_stl --update-damping`, then copy
the resulting MAT-file here.

The source geometry is credited to the OpenFOAM Foundation, via Ddcampayo's
[Wikimedia Commons upload](https://commons.wikimedia.org/wiki/File:Wigley_hull.stl).
The adapted mesh is distributed in MSS-Capytaine under CC BY-SA 4.0; see its
[source and licence notice](https://github.com/cybergalactic/MSS-Capytaine/blob/main/vessels_capytaine/wigley_stl/README.md#source-and-licence).
Geometry files are not included in this MSS directory.
