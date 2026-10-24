# Capytaine vessel data

This directory contains pre-generated hydrodynamic vessel structures exported
by [MSS-Capytaine](https://github.com/cybergalactic/MSS-Capytaine). The files
use the common MSS `vessel` format and can be loaded directly in MATLAB or GNU
Octave; Python and Capytaine are needed only to regenerate them.

| Case | Type | File |
| --- | --- | --- |
| [`testShip`](testShip/Test-Ship.md) | Synthetic surface monohull | `testShip/testShip.mat` |
| [`LAUV_marie`](LAUV_marie/LAUV-Marie.md) | Idealized submerged LAUV Marie model | `LAUV_marie/LAUV_marie.mat` |

The data use MSS forward-starboard-down body axes and are referenced to the
center of gravity. Both cases are demonstrations of the Capytaine-to-MSS
workflow rather than validated vessel designs.

LAUV Marie also has a dedicated standalone craft model,
[`lauvMarie.m`](../../CRAFT/AUV/lauvMarie.m), and wave-enabled simulation,
[`SIMlauvMarie.m`](../../CRAFT/AUV/SIMlauvMarie.m). These add LAUV-specific
nonlinear and calibratable terms, controls, guidance, and depth attenuation of
wave forces. Use `SIMlauvMarie` for the dedicated LAUV model. Its MAT-file
can also be loaded manually in `SIMhydroVessel` for a generic hydrodynamic
test. The GUI initially loads `LAUV_marie.mat` as a generic custom model.

The generated files also contain the linear viscous-damping inputs selected
in each MSS-Capytaine `config.json` under `vessel.powerBased`. For these two
cases, `hydroVesselConfig.m` uses the stored values as the simulation defaults.
Surface vessels use `T_126 = [T_surge T_sway T_yaw]` in seconds and
`delta_zeta_345`; submerged vehicles use `T_1236` in seconds and
`delta_zeta_45`. Surface time constants replace the former `kappa_126`
percentages and must not require negative additional viscous damping.
Commented assignments in the corresponding vessel cases can be enabled to
override them. WAMIT and ShipX files do not currently store these inputs, so
their active defaults remain in `hydroVesselConfig.m`.

## Regenerating the data

Clone MSS-Capytaine separately, install its Python dependencies, and run the
desired case from the MSS-Capytaine repository root:

```sh
python main.py testShip
python main.py LAUV_marie
```

MSS-Capytaine writes each generated vessel file under
`vessels_capytaine/capytaine_offsets_csv/<case>/results/` for these two cases. Copy the resulting `.mat` file into the
matching directory here to refresh the pre-generated MSS catalogue.
Only the vessel MAT-file is required. Verify geometry in MSS-Capytaine before
exporting; the MSS test GUI does not require geometry files or plot the mesh.

For a newly developed vessel, run `SIMhydroVessel` and click
**Load vessel MAT-file** to open the generated vessel file directly from
`MSS-Capytaine/vessels_capytaine/<family>/<case>/results/`.
The GUI retains TestShip as a demonstration and imported ShipX/WAMIT data, labeled by solver. The GUI loads `vessel.powerBased` and presents
conservative simulation and controller defaults for review and editing. In
particular, the generic constant surge force is zero and must be selected by
the user when forward propulsion is required. Configure the remaining settings,
then click **Run simulation**. Copying the custom file into the MSS catalogue
is optional.
