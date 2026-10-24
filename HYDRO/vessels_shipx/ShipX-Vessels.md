# ShipX vessel data

This directory contains two ShipX (VERES) data sets:

- `s175`: S-175 container ship
- `supply`: offshore supply vessel

Each vessel is processed in two steps. First, `veres2vessel` reads the
ShipX files and creates the frequency-dependent MSS vessel structure.
Second, `computeManeuveringModel` forms the constant power-based maneuvering
model for a selected sea state.

Run the commands from the directory containing the ShipX files for the vessel
being processed. The ShipX source files in both directories use the basename
`input`. When `veres2vessel` asks for the vessel name used for the output file,
enter `s175` or `supply` as indicated below.

## S-175 container ship

From the `s175` directory, run:

```matlab
vessel = veres2vessel('input','1111');
```

At the output-name prompt, enter:

```text
s175
```

Then compute the power-based model:

```matlab
vessel = computeManeuveringModel(vessel,omega_p,T_126,delta_zeta_345,1);
```

## Offshore supply vessel

From the `supply` directory, run:

```matlab
vessel = veres2vessel('input','1111');
```

At the output-name prompt, enter:

```text
supply
```

Then compute the power-based model:

```matlab
vessel = computeManeuveringModel(vessel,omega_p,T_126,delta_zeta_345,1);
```

## Plot commands

After running the two processing commands for a vessel, the following commands
can be used to inspect the results:

```matlab
% Frequency-dependent hydrodynamic matrices
plotABC(vessel,'A');
plotABC(vessel,'B');

% Potential and power-based damping matrices
plotBv(vessel);

% Motion and force response amplitude operators
plotTF(vessel,'motion','rads',1);
plotTF(vessel,'force','rads',1);

% Second-order wave-drift forces
plotWD(vessel,'rads',1);
```

The final argument `1` selects the zero-speed ShipX data set. The power-based
damping plot requires `computeManeuveringModel` to have been run first.

ShipX also supplies Ikeda roll damping. This remains available separately as
`vessel.roll.Bv44`; it is not the power-based viscous damping matrix stored as
`vessel.powerBased.Bv` by `computeManeuveringModel`.

## Power-based model arguments

The second processing command is:

```matlab
vessel = computeManeuveringModel(vessel,omega_p,T_126,delta_zeta_345,plotFlag);
```

- `vessel` is the frequency-dependent structure returned by `veres2vessel`.
- `omega_p` is the wave-spectrum peak frequency in rad/s.
- `T_126` contains positive total velocity-decay time constants in seconds
  for surge, sway, and yaw (DOFs 1, 2, and 6). There is no universal default.
- `delta_zeta_345` contains the additional viscous damping ratios for heave,
  roll, and pitch (DOFs 3, 4, and 5). The default is `[0 0.1 0]`.
- `plotFlag` is `1` to plot the frequency-dependent matrices together with
  their power-based equivalent values, or `0` for no plots.

The surface-vessel presets in `hydroVesselConfig.m` are feasible starting
values at `omega_p = 1.0` rad/s and zero speed:

| Vessel | `T_126` (s) | `delta_zeta_345` |
| --- | --- | --- |
| S-175 | `[11.9 2.0 1.8]` | `[0 0.2 0]` |
| Supply | `[6.4 2.5 3.3]` | `[0 0.1 0]` |

For example, for the supply vessel:

```matlab
omega_p = 1.0;
T_126 = [6.4 2.5 3.3];
delta_zeta_345 = [0 0.1 0];

vessel = computeManeuveringModel(vessel,omega_p,T_126,delta_zeta_345,1);
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

The result contains the constant matrices and selected inputs under
`vessel.powerBased`, together with the maneuvering-model matrices `M`, `D`, and
`G`. Changing `omega_p`, `T_126`, or `delta_zeta_345` only requires rerunning
`computeManeuveringModel`; the ShipX import does not need to be repeated.

## ShipX plot flag

The optional string passed to `veres2vessel` controls its plots:

- `'1000'`: added-mass and potential-damping matrices
- `'0100'`: force RAOs
- `'0010'`: motion RAOs
- `'0001'`: wave-drift forces
- `'0000'`: no plots
- `'1111'`: all plots

If omitted, the default is `'1000'`.

`veres2vessel` reads `*.re1`, `*.re2`, `*.re7`, `*.re8`, and `*.hyd`, and saves
the imported frequency-dependent vessel structure using the name entered at the
prompt. The subsequent call to `computeManeuveringModel` updates the structure in
the MATLAB or GNU Octave workspace. Save it again if the derived power-based
model should be retained on disk, for example:

```matlab
save('s175.mat','vessel');
% or
save('supply.mat','vessel');
```
