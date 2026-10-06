# Capytaine vessel data

This directory contains pre-generated hydrodynamic vessel structures exported
by [MSS-Capytaine](https://github.com/cybergalactic/MSS-Capytaine). The files
use the common MSS `vessel` format and can be loaded directly in MATLAB or GNU
Octave; Python and Capytaine are needed only to regenerate them.

| Case | Type | File |
| --- | --- | --- |
| [`testShip`](testShip/) | Synthetic surface monohull | `testShip/testShip.mat` |
| [`LAUV_marie`](LAUV_marie/) | Idealized submerged LAUV Marie model | `LAUV_marie/LAUV_marie.mat` |

The data use MSS forward-starboard-down body axes and are referenced to the
center of gravity. Both cases are demonstrations of the Capytaine-to-MSS
workflow rather than validated vessel designs.

The generated files also contain the linear viscous-damping inputs selected
in each MSS-Capytaine `config.json` under `vessel.powerBased`. For these two
cases, `hydroVesselConfig.m` uses the stored values as the simulation defaults.
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
`vessels_capytaine/<case>/results/`. Copy the resulting `.mat` file into the
matching directory here to refresh the pre-generated MSS catalogue.
