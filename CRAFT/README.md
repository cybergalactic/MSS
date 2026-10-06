# CRAFT Library

CRAFT contains ready-to-run mathematical models and editable time-domain simulation scripts for ships, autonomous underwater vehicles (AUVs), uncrewed surface vehicles (USVs), and floating structures. It also provides a unified workflow for simulating user-defined 6-DOF marine craft from ShipX, WAMIT, or Capytaine seakeeping data. The files are compatible with MATLAB and GNU Octave.

## Modeling Workflows

### Established craft models

Each established model is paired with an editable `SIM*.m` script that demonstrates guidance, control, maneuvering, or dynamic positioning.

| Type | Model | Simulation |
| --- | --- | --- |
| AUV | [`DSRV.m`](AUV/models/DSRV.m) | [`SIMdsrv.m`](AUV/SIMdsrv.m) |
| AUV | [`npsauv.m`](AUV/models/npsauv.m) | [`SIMnpsauv.m`](AUV/SIMnpsauv.m) |
| AUV | [`remus100.m`](AUV/models/remus100.m) | [`SIMremus100.m`](AUV/SIMremus100.m) |
| USV | [`otter.m`](USV/models/otter.m) | [`SIMotter.m`](USV/SIMotter.m) |
| Ship | [`clarke83.m`](SHIP/models/clarke83.m) | [`SIMclarke83.m`](SHIP/SIMclarke83.m) |
| Ship | [`container.m`](SHIP/models/container.m) and [`Lcontainer.m`](SHIP/models/Lcontainer.m) | [`SIMcontainer.m`](SHIP/SIMcontainer.m) |
| Ship | [`frigate.m`](SHIP/models/frigate.m) | [`SIMfrigate.m`](SHIP/SIMfrigate.m) |
| Ship | [`mariner.m`](SHIP/models/mariner.m) | [`SIMmariner.m`](SHIP/SIMmariner.m) |
| Ship | [`navalvessel.m`](SHIP/models/navalvessel.m) | [`SIMnavalvessel.m`](SHIP/SIMnavalvessel.m) |
| Ship | [`osv.m`](SHIP/models/osv.m) | [`SIMosv.m`](SHIP/SIMosv.m) |
| Floating structure | [`rig.m`](SHIP/models/semisubModels/rig.m) | [`SIMsemisub.m`](SHIP/SIMsemisub.m) |
| Ship | [`supply.m`](SHIP/models/supply.m) | [`SIMsupply.m`](SHIP/SIMsupply.m) |
| Ship | [`tanker.m`](SHIP/models/tanker.m) | [`SIMtanker.m`](SHIP/SIMtanker.m) |
| Craft | [`zeefakkel.m`](SHIP/models/zeefakkel.m) | [`SIMzeefakkel.m`](SHIP/SIMzeefakkel.m) |

### User-defined hydrodynamic craft

[`hydroVessel.m`](hydroVessel.m) implements a 12-state, nonlinear 6-DOF model using the common MSS `vessel` structure. [`SIMhydroVessel.m`](SIMhydroVessel.m) is an editable simulation template that can be adapted to a user-defined USV, AUV, ship, or floating structure.

The simulator includes a basic GUI for choosing a vessel and configuring irregular seas, wave spectra, directional spreading, first-order wave forces from force RAOs, ocean currents, viscous damping, initial conditions, and either heading-autopilot or dynamic-positioning control.

Hydrodynamic vessel data and processing instructions are organized by source:

- [ShipX vessel data](../HYDRO/vessels_shipx/)
- [WAMIT vessel data](../HYDRO/vessels_wamit/)
- [Capytaine vessel data](../HYDRO/vessels_capytaine/)

ShipX and WAMIT are commercial programs. The open-source [MSS-Capytaine](https://github.com/cybergalactic/MSS-Capytaine) add-on uses Capytaine to generate compatible vessel structures without commercial seakeeping software.

## Getting Started

From the MSS root directory, update the MATLAB path and run a simulation:

```matlab
mssPath
SIMremus100
```

MATLAB users can also display the concise command-window catalogue:

```matlab
help CRAFT
```

To use the hydrodynamic workflow, run:

```matlab
SIMhydroVessel
```

The configuration and GUI helpers are located in [`utils_hydroVessel`](utils_hydroVessel/).
