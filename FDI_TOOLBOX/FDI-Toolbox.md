# Frequency-Domain Identification Toolbox

For marine craft, the FDI toolbox uses frequency-dependent added mass
`A(ω)` and radiation damping `B(ω)` to fit retardation functions that
represent fluid memory in the Cummins equation. It can also estimate the
infinite-frequency added mass `A∞` when it is unavailable. Fitting and
checking these models takes more effort, so the default MSS maneuvering
workflow uses the power-based method of [Fossen (2025)](https://doi.org/10.1016/j.oceaneng.2025.122494).
The function [`computeManeuveringModel.m`](../HYDRO/computeManeuveringModel.m)
averages `A(ω)` and `B(ω)` over a reference wave spectrum to obtain
constant equivalent matrices. Use FDI when the application calls for an
explicit fluid-memory model, as in published studies of a three-tether
wave energy converter and the wave-induced bending of a multibody subsea
tanker. See [Example applications](#example-applications).

## Demonstrations

- [`Demo_FDIRadMod_WA.m`](Demo_FDIRadMod_WA.m): identification when the infinite-frequency added mass is available.
- [`Demo_FDIRadMod_NA.m`](Demo_FDIRadMod_NA.m): identification when the infinite-frequency added mass must be estimated.

Both demonstrations load the shared [WAMIT FPSO dataset](../HYDRO/vessels_wamit/fpso.mat) by name from the MSS path. FDI does not keep a separate vessel copy; its identification inputs are the same added-mass and radiation-damping data used previously.

## Main Functions

| Function | Purpose |
| --- | --- |
| `FDIRadMod.m` | Identify a stable SISO radiation-force transfer function for a selected hydrodynamic coupling. |
| `EditAB.m` | Interactively select a frequency range and remove outliers from added-mass and damping data. |
| `ident_retardation_FD.m` | Fit a parametric radiation-retardation model from complex frequency-response data. |
| `ident_retardation_FDna.m` | Identify a parametric model using the alternative added-mass formulation. |
| `fit_siso_fresp.m` | Fit a continuous SISO transfer function to frequency-response data. |

## Getting Started

Refresh the MSS search path, then run one of the demonstrations from the
MATLAB Command Window:

```matlab
mssPath
Demo_FDIRadMod_WA
```

If your saved path includes the former `FDI/` folder, run `mssPath` again in
MATLAB. In GNU Octave, refresh the saved path using the
[Octave installation guide](../How%20to%20install%20MSS%20for%20GNU%20OCTAVE.md).

## Tutorial and papers

- [Identification of Seakeeping Models from Frequency Response Data](DOCS/FDI%20Tutorial.pdf).
- T. Perez and T. I. Fossen (2008). [Time- vs Frequency-Domain Identification of Parametric Radiation Force Models for Marine Structures at Zero Speed](DOCS/PerezFossen%20MIC-29%281%29%202008.pdf). *Modeling, Identification and Control*, 29(1), 1–19. <https://doi.org/10.4173/mic.2008.1.1>
- T. Perez and T. I. Fossen (2008). [Joint Identification of Infinite-Frequency Added Mass and Fluid-Memory Models of Marine Structures](DOCS/PerezFossen%20MIC-29%283%29%202008.pdf). *Modeling, Identification and Control*, 29(3), 93–102. <https://doi.org/10.4173/mic.2008.3.2>

T. Perez and T. I. Fossen (2009). “A MATLAB Tool for Parametric Identification of Radiation-Force Models of Ships and Offshore Structures.” *Modeling, Identification and Control*, 30(1), 1–15. <https://doi.org/10.4173/mic.2009.1.1>

## Example applications

- N. Y. Sergiienko (2018). *Three-Tether Wave Energy Converter: Hydrodynamic
  Modelling, Performance Assessment and Control*. PhD thesis, University of
  Adelaide, §3.4.2. Uses MSS FDI to model the converter's radiation forces
  in the time domain. [Thesis PDF](https://digital.library.adelaide.edu.au/server/api/core/bitstreams/6b2f442d-6458-4df2-b14b-400d6dc1c981/content).
- K. S. Patel (2023). *Bending Moment and Efficient Fatigue Assessment in a
  Subsea Shuttle Tanker Under the Effect of Waves*. Master's thesis,
  University of Stavanger, §3.3 and Appendix B. Fits radiation models for a
  seven-body representation of the tanker. [Thesis PDF](https://uis.brage.unit.no/uis-xmlui/bitstream/handle/11250/3084649/no.uis%3Ainspera%3A137508920%3A90468133.pdf?isAllowed=y&sequence=1).
