# Capytaine vessel data

This directory contains pre-generated hydrodynamic vessel structures exported
by [MSS-Capytaine](https://github.com/cybergalactic/MSS-Capytaine). The files
use the common MSS `vessel` format and can be loaded directly in MATLAB or GNU
Octave; Python and Capytaine are needed only to regenerate them.

| Case | Type | File |
| --- | --- | --- |
| [`testShip`](testShip/Test-Ship.md) | Synthetic surface monohull | `testShip/testShip.mat` |
| [`LAUV_marie`](LAUV_marie/LAUV-Marie.md) | Idealized submerged LAUV Marie model | `LAUV_marie/LAUV_marie.mat` |
| [Supply vessel](supply/Supply.md) | Surface monohull | [`supply/supply.mat`](supply/supply.mat) |
| [Wigley](wigley/Wigley.md) | Benchmark hull | [`wigley/wigley.mat`](wigley/wigley.mat) |
| [S175](S175/S175.md) | Container ship | [`S175/s175.mat`](S175/s175.mat) |
| [Otter USV](#otter-usv) | Twin-hull surface craft, 55 kg and 80 kg load conditions | [`otter/otter_55kg.mat`](otter/otter_55kg.mat), [`otter/otter_80kg.mat`](otter/otter_80kg.mat) |

The data use MSS forward-starboard-down body axes and are referenced to the
center of gravity. The cases demonstrate the Capytaine-to-MSS workflow;
their hydrodynamic data and dedicated craft models require validation against
measurements for quantitative vessel-response prediction.

The included vessel MAT-files contain signed mean-drift coefficients in
`vessel.drift.w` and `vessel.drift.coeff{1:3}`, and complete damping inputs
in `vessel.powerBased`.
The copies match the source files byte for byte. Refresh them after updating
the source catalogue.

LAUV Marie also has a dedicated standalone craft model,
[`lauvMarie.m`](../../CRAFT/AUV/models/lauvMarie.m), and wave-enabled simulation,
[`SIMlauvMarie.m`](../../CRAFT/AUV/SIMlauvMarie.m). These add LAUV-specific
nonlinear and calibratable terms, controls, guidance, and depth attenuation of
wave forces. Use `SIMlauvMarie` for the dedicated LAUV model. Its MAT-file
can also be loaded manually in `SIMhydroVessel` for a generic hydrodynamic
test. The GUI initially loads `LAUV_marie.mat` as a generic custom model.

The generated files also contain the linear viscous-damping inputs selected
in each MSS-Capytaine `config.json` under `vessel.powerBased`.
`hydroVesselConfig.m` uses these stored values as the simulation defaults.
Surface vessels use `T_126 = [T_surge T_sway T_yaw]` in seconds and
`delta_zeta_345`; submerged vehicles use `T_1236` in seconds and
`delta_zeta_45`. Surface time constants replace the former `kappa_126`
percentages and must not require negative additional viscous damping.
Simulation-specific overrides can be set in `hydroVesselConfig.m`. All three
WAMIT vessels also store their damping inputs, collected by `wamit2vessel`. Power-based averaging defaults to a JONSWAP
reference with peak 1.0 rad/s and gamma 3.3. An explicitly supplied peak
overrides this reference; the excitation spectrum is configured separately.

These Capytaine files contain a single zero-speed hydrodynamic dataset.
They do not replace speed-dependent radiation matrices for forward-speed
studies. The 2025 comparison retains its legacy S175 file separately under
[`mssExamples/data`](../../mssExamples/data/README.md).

## Otter USV

The Otter files contain frequency-dependent added mass and radiation damping,
hydrostatic restoring, and force and motion RAOs for the 55 kg hull and the
nominal 80 kg loaded condition. The dedicated
[`otter.m`](../../CRAFT/USV/models/otter.m) model uses `otter_80kg.mat`, with
power-based equivalent matrices computed for a fixed reference spectrum.
It transforms the matrices from the dataset CG to the body-frame origin CO
and adds the vessel-specific nonlinear maneuvering terms.

[`SIMotter.m`](../../CRAFT/USV/SIMotter.m) provides heading/path-following
control, currents, and optional directional first-order wave loads. The same
80 kg force RAOs can be applied to the new model or to the first-principles
[`otter_legacy.m`](../../CRAFT/USV/models/otter_legacy.m) using the GUI's
legacy-model checkbox. Wave loads are transformed from the dataset CG to CO
before they enter the dynamics. Changing the payload does not regenerate the
nominal hydrodynamic data or hydrostatic stiffness.

## Regenerating the data

Clone MSS-Capytaine separately, install its Python dependencies, and run the
desired case from the MSS-Capytaine repository root:

```sh
python main.py testShip_csv
python main.py LAUV_marie_csv
python main.py supply_mgf
python main.py s175_mgf
```

MSS-Capytaine 1.2 writes each main vessel MAT-file beside its `config.json`
under `vessels_capytaine/<vessel>_<format>/`; supporting files remain in
`results/`. Copy the resulting `.mat` file into the
matching directory here to refresh the pre-generated MSS catalogue.
Only the vessel MAT-file is required. Verify geometry in MSS-Capytaine before
exporting; the MSS test GUI does not require geometry files or plot the mesh.

For a newly developed vessel, run `SIMhydroVessel` and click
**Load vessel MAT-file** to open the generated vessel file directly from
`MSS-Capytaine/vessels_capytaine/<vessel>_<format>/`.
The GUI groups Test ship, Supply vessel and S175 under Capytaine, and Tanker,
FPSO and Semisubmersible under WAMIT. The GUI loads `vessel.powerBased` and presents
conservative simulation and controller defaults for review and editing. In
particular, the generic constant surge force is zero and must be selected by
the user when forward propulsion is required. Configure the remaining settings,
then click **Run simulation**. Copying the custom file into the MSS catalogue
is optional.

To change only the saved viscous-damping inputs after editing a case JSON,
run `python main.py supply_mgf --update-damping` (substitute the desired case).
This updates the existing MAT-file without rerunning Capytaine. Copy the
updated file into MSS afterwards; equivalent matrices are computed when MSS
loads and processes the vessel.

MSS includes the WAMIT tanker as `tanker.mat`. The Capytaine tanker remains
in MSS-Capytaine as a validation case against WAMIT and is not copied into
the MSS catalogue.
