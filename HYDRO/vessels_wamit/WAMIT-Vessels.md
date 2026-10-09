# WAMIT vessel data

This directory contains three WAMIT data sets:

- `tanker`
- `semisub`
- `fpso`

Each vessel is processed in two steps. First, `wamit2vessel` reads the WAMIT files and creates the frequency-dependent MSS vessel structure.
Second, `computeManeuveringModel` forms the constant power-based maneuvering model for a selected sea state.

Run the commands from the directory containing the WAMIT files for the vessel being processed.

## Tanker

```matlab
vessel = wamit2vessel('tanker',10,246,46);
vessel = computeManeuveringModel(vessel,omega_p,kappa_126,delta_zeta_345,1);
```

To display all WAMIT plots during import, use:

```matlab
vessel = wamit2vessel('tanker',10,246,46,'1111');
vessel = computeManeuveringModel(vessel,omega_p,kappa_126,delta_zeta_345,1);
```

## Semisubmersible

```matlab
vessel = wamit2vessel('semisub',21,115,80,'1111');
vessel = computeManeuveringModel(vessel,omega_p,kappa_126,delta_zeta_345,1);
```

## FPSO

```matlab
vessel = wamit2vessel('fpso',12,200,44,'1111');
vessel = computeManeuveringModel(vessel,omega_p,kappa_126,delta_zeta_345,1);
```

## Plot commands

After running the two processing commands for a vessel, the following commands can be used to inspect the results:

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

The final argument `1` selects the zero-speed WAMIT data set. The
power-based damping plot requires `computeManeuveringModel` to have been
run first.

## Power-based model arguments

The second processing command is:

```matlab
vessel = computeManeuveringModel(vessel,omega_p,kappa_126,delta_zeta_345,plotFlag);
```

- `vessel` is the frequency-dependent structure returned by
  `wamit2vessel`.
- `omega_p` is the wave-spectrum peak frequency in rad/s.
- `kappa_126` contains the relative viscous-damping increments for surge, sway, and yaw (DOFs 1, 2, and 6). The default is `[0.05 0.05 0.05]`.
- `delta_zeta_345` contains the additional viscous damping ratios for heave, roll, and pitch (DOFs 3, 4, and 5). The default is `[0 0.1 0]`.
- `plotFlag` is `1` to plot the frequency-dependent matrices together with their power-based equivalent values, or `0` for no plots.

For example:

```matlab
omega_p = 1.0;
kappa_126 = [0.05 0.05 0.05];
delta_zeta_345 = [0 0.1 0];

vessel = computeManeuveringModel(vessel,omega_p,kappa_126,delta_zeta_345,1);
```

To enter the damping values interactively and accept the displayed defaults by pressing Return, use:

```matlab
vessel = computeManeuveringModel(vessel,omega_p,[],[],1);
```

The result contains the constant matrices and selected inputs under `vessel.powerBased`, together with the maneuvering-model matrices `M`,
`D`, and `G`. Changing `omega_p`, `kappa_126`, or `delta_zeta_345` only requires rerunning `computeManeuveringModel`; the WAMIT import does not
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

`wamit2vessel` saves the imported frequency-dependent vessel data. The subsequent call to `computeManeuveringModel` updates the vessel structure
in the MATLAB or GNU Octave workspace. Save it again if the derived power-based model should be retained on disk.
