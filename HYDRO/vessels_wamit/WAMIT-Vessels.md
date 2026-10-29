# WAMIT vessel data

This catalogue has three ready-to-load MAT-files beside `wamit2vessel.m`:

- `tanker.mat`
- `semisub.mat`
- `fpso.mat`

Each vessel is processed in two steps. First, `wamit2vessel` reads the WAMIT files, collects editable damping inputs and saves the frequency-dependent MSS vessel structure with `vessel.powerBased` metadata.
Second, `computeManeuveringModel` forms the constant power-based maneuvering model for a selected sea state.

The raw WAMIT files remain in the `tanker_wamitData`,
`semisub_wamitData`, and `fpso_wamitData` subfolders.
The importer always saves the MAT-file beside itself, even when called from
another working directory. The three MAT-files are then available by name on
the MSS path, for example `load tanker.mat`.

## Regenerate the three MAT-files

Run these calls from `HYDRO/vessels_wamit`. They use the stored damping defaults
and do not open plots or request input:

```matlab
wamit2vessel(fullfile('tanker_wamitData','tanker'),10,246,46,'0000', ...
    [49.6 3.6 4.5],[0 0.1 0],[5 0 0]);
wamit2vessel(fullfile('semisub_wamitData','semisub'),21,115,80,'0000', ...
    [6.1 3.1 3.1],[0 0.1 0],[5 0 0]);
wamit2vessel(fullfile('fpso_wamitData','fpso'),12,200,44,'0000', ...
    [22.7 3.1 3.4],[0 0.1 0],[5 0 0]);
```

For interactive damping choices, omit the last three vector arguments; Return
accepts each proposal. To display plots during import, change `'0000'` to
`'1111'`.

## Current vessel structure

The importer and bundled MAT-files use `vessel.drift.w` and
`vessel.drift.coeff{1:3}` for signed surge, sway and yaw mean-drift
coefficients. These are loads per squared incident-wave amplitude, in N/m²
and Nm/m². They are read from the signed real part of the WAMIT `.8` or `.9`
diagonal output and converted to MSS forward-starboard-down axes. The
180–360 degree half-plane uses port-starboard symmetry.
Scaling follows the [WAMIT mean-drift normalization](https://www.wamit.com/manual7.x/html/wamit_v75manualch3.html):
multiply nondimensional forces by `rho*g*ULEN` and moments by `rho*g*ULEN^2`
to obtain coefficients per squared wave amplitude.

The WAMIT body origin is retained, recorded in `vessel.drift.reference`;
these moments have not been translated to the center of gravity. Capytaine
exports use their documented CG reference. The saved RAO and drift arrays contain
one speed case, consistent with `vessel.velocities = 0`.

The former `vessel.driftfrc.amp` field has been removed from the bundled
files. `plotWD` uses only the current `vessel.drift` structure. The migration
does not change rigid-body properties, radiation matrices, first-order RAO
values or damping inputs, and requires no new WAMIT solve. Existing `*ABC.mat`
files are separate radiation models and contain no drift fields.

## Plot commands

For example, load a saved vessel and compute the equivalent model before
inspecting the results:

```matlab
load tanker.mat
vessel = computeManeuveringModel(vessel);

% Frequency-dependent hydrodynamic matrices
plotABC(vessel,'A');
plotABC(vessel,'B');

% Potential and power-based damping matrices
plotBv(vessel);

% Motion and force response amplitude operators
plotRAO(vessel);

% Second-order wave-drift forces
plotWD(vessel,'rads',1);
```

`plotRAO` shows both motion and force RAOs for the default headings
0, 30, 60 and 90 degrees, using the first (zero-speed) data set.
The final argument `1` in `plotWD` also selects zero speed. The power-based
damping plot requires `computeManeuveringModel` to have been run first.
Use `plotTF` only when you need separate motion or force plots, a different
frequency unit, or another speed case.

## Power-based model arguments

The second processing command is:

```matlab
vessel = computeManeuveringModel(vessel,omega_p,T_126,delta_zeta_345,plotFlag);
```

- `vessel` is the frequency-dependent structure returned by
  `wamit2vessel`.
- `omega_p` is the reference JONSWAP peak frequency in rad/s; omitted or empty uses 1.0 rad/s, with gamma 3.3.
- `T_126` contains positive total velocity-decay time constants in seconds for surge, sway, and yaw (DOFs 1, 2, and 6). There is no universal default.
- `delta_zeta_345` contains the additional viscous damping ratios for heave, roll, and pitch (DOFs 3, 4, and 5). The default is `[0 0.1 0]`.
- `plotFlag` is `1` to plot the frequency-dependent matrices together with their power-based equivalent values, or `0` for no plots.

All three bundled WAMIT vessels store damping inputs in `vessel.powerBased`.
`wamit2vessel` proposes the following defaults at the JONSWAP reference
`omega_p = 1.0` rad/s and zero speed. Press Return to accept each proposal,
or enter another vector:

| Vessel | `T_126` (s) |
| --- | --- |
| Tanker | `[49.6 3.6 4.5]` |
| Semisubmersible | `[6.1 3.1 3.1]` |
| FPSO | `[22.7 3.1 3.4]` |

The proposed damping-ratio increment is `delta_zeta_345 = [0 0.1 0]`,
and the nonlinear angular multipliers are `nonlinear_456 = [5 0 0]`.
The importer checks linear damping feasibility before saving. New vessel
names require explicit time constants because there is no universal default.

For a scripted import without prompts, supply all three vectors:

```matlab
vessel = wamit2vessel('tanker',10,246,46,'0000', ...
    [49.6 3.6 4.5],[0 0.1 0],[5 0 0]);
```

The saved `powerBased` field contains these inputs and the reference spectrum;
`A_eq`, `B_eq`, `Bv`, `M`, `D` and `G` are calculated separately at runtime.

For example, for the tanker:

```matlab
omega_p = 1.0;
T_126 = [49.6 3.6 4.5];
delta_zeta_345 = [0 0.1 0];

vessel = computeManeuveringModel(vessel,[],vessel.powerBased.T_126,vessel.powerBased.delta_zeta_345,1);
```

To enter damping values interactively, use the following call. Enter all three
`T_126` values; Return accepts the default only for `delta_zeta_345`:

```matlab
vessel = computeManeuveringModel(vessel,omega_p,[],[],1);
```

For the diagonal equivalent model, `Bv(i,i) = M(i,i)/T_i - B_eq(i,i)`
for DOFs 1, 2, and 6. Smaller time constants increase damping. Targets that
require negative additional viscous damping are rejected; for positive
`B_eq(i,i)`, choose `T_i <= M(i,i)/B_eq(i,i)`. Coupled decay times can differ.
Recheck feasibility when changing the sea state or speed.

The result contains the constant matrices and selected inputs under `vessel.powerBased`, together with the maneuvering-model matrices `M`,
`D`, and `G`. Changing `omega_p`, `T_126`, or `delta_zeta_345` only requires rerunning `computeManeuveringModel`; the WAMIT import does not
need to be repeated.

## WAMIT plot flag

The optional string passed to `wamit2vessel` controls its plots:

- `'1000'`: added-mass and potential-damping matrices
- `'0100'`: force RAOs
- `'0010'`: motion RAOs
- `'0001'`: wave-drift forces
- `'0000'`: no plots
- `'1111'`: all plots

If omitted, the default is `'1000'`.

`wamit2vessel` saves the imported frequency-dependent vessel data and chosen damping inputs. The subsequent call to `computeManeuveringModel` updates the vessel structure
in the MATLAB or GNU Octave workspace. Save it again if the derived power-based model should be retained on disk.

MSS includes only the WAMIT tanker dataset, named `tanker.mat`. Load it with
`load tanker`; `load(resolveHydroVesselFile('tanker.mat'),'vessel')` also
selects this dataset. The Capytaine tanker is kept in MSS-Capytaine for
validation against the WAMIT results.
