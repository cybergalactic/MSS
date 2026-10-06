# Frequency-Domain Identification Toolbox

The FDI toolbox identifies parametric radiation-force and fluid-memory models of ships, offshore structures, and wave-energy systems from frequency-dependent added-mass and potential-damping data.

## Demonstrations

- [`Demo_FDIRadMod_WA.m`](Demo_FDIRadMod_WA.m): identification when the infinite-frequency added mass is available.
- [`Demo_FDIRadMod_NA.m`](Demo_FDIRadMod_NA.m): identification when the infinite-frequency added mass must be estimated.

Both demonstrations use the included `fpso.mat` vessel data.

## Main Functions

| Function | Purpose |
| --- | --- |
| `FDIRadMod.m` | Identify a stable SISO radiation-force transfer function for a selected hydrodynamic coupling. |
| `EditAB.m` | Interactively select a frequency range and remove outliers from added-mass and damping data. |
| `ident_retardation_FD.m` | Fit a parametric radiation-retardation model from complex frequency-response data. |
| `ident_retardation_FDna.m` | Identify a parametric model using the alternative added-mass formulation. |
| `fit_siso_fresp.m` | Fit a continuous SISO transfer function to frequency-response data. |

## Getting Started

Add MSS and its subfolders to the MATLAB path, change to the `FDI` directory, and run one of the demonstrations:

```matlab
mssPath
Demo_FDIRadMod_WA
```

The detailed tutorial and supporting publications are available in the [`documentation/FDI Identification of seakeeping models from frequency response data`](../documentation/FDI%20Identification%20of%20seakeeping%20models%20from%20frequency%20response%20data/) directory.

## Reference

T. Perez and T. I. Fossen (2009). “A MATLAB Tool for Parametric Identification of Radiation-Force Models of Ships and Offshore Structures.” *Modeling, Identification and Control*, 30(1), 1–15. <https://doi.org/10.4173/mic.2009.1.1>
