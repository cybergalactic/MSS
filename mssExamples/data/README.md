# Example-specific vessel data

`s175_speedDependent.mat` retains the legacy ShipX S175 hydrodynamics used
by [`exManeuveringModel`](../exManeuveringModel.m) and the 2025 power-based
averaging study. It contains five speed cases: 0, 2.5722, 5.1444, 7.7167 and
10.2889 m/s. The example selects the third speed case, 5.1444 m/s.

This dataset is separate from the zero-speed Capytaine `s175.mat` in
`HYDRO/vessels_capytaine/S175/` and is omitted from the `mssCraft`
hydrodynamic catalogue. The ShipX importer and general catalogue have been
removed; this retained file needs neither ShipX nor its importer.

The comparison uses the same JONSWAP reference at 1.0 rad/s, gamma 3.3,
and additive viscous damping in both methods. Its saved damping inputs are
`T_126 = [11.4 1.8 1.6]` s and `delta_zeta_345 = [0 0.1 0]`.
The default run is 300 s so the response can be compared after startup
transients have decayed. Radiation memory is retained for 50 s.

The retained file also uses the current `vessel.drift.w` /
`vessel.drift.coeff{1:3}` layout. Empty speed padding has been removed from
RAO and drift arrays; the five physical speed cases are preserved. Opposite
headings follow surge-even, sway-odd and yaw-odd symmetry. Radiation matrices,
first-order RAO values and saved damping inputs are unchanged.
