# Getting Started with MSS

Install the Marine Systems Simulator (MSS), run a ship heading-control
simulation, and find the models and examples for your application.

## 1. Install MSS

- **MATLAB.** Follow the [MATLAB installation guide](../How%20to%20install%20MSS%20for%20MATLAB.md).
  After installing MSS or adding its root folder to the MATLAB path, run
  `mssPath` to add the subfolders and save the path.
- **GNU Octave.** Follow the [GNU Octave installation guide](../How%20to%20install%20MSS%20for%20GNU%20OCTAVE.md)
  to start Octave with Qt graphics and add the MSS folder and its subfolders
  to the search path.

In the Command Window, check that both commands resolve to your MSS
installation:

```matlab
which mssCraft
which SIMfrigate
```

If either command is missing or points to a different copy, revisit the path
setup in the installation guide for your environment.

## 2. Run a first simulation

Enter:

```matlab
SIMfrigate
```

The [frigate simulation](../CRAFT/SHIP/SIMfrigate.m) runs with preset
parameters and opens plots of heading, yaw rate, and rudder angle. Compare
the simulated heading with the desired heading. The heading command starts
at 30° and changes to 70° just after 100 s.

To try path following next, run `SIMotter`, select **ALOS path following:
straight lines and waypoint switching**, and click **Start simulation**.
The [Otter simulation](../CRAFT/USV/SIMotter.m) lets you change the guidance
method, current, waves, and simulation time in the setup window.

## 3. Find the right model or example

Run `mssCraft` to list the bundled hydrodynamic MAT-files and the CRAFT
dynamics models with their companion simulations. Enter a simulation name
without the `.m` extension to run it.

| Goal | Where to start |
| --- | --- |
| Simulate a ship, AUV, USV, or floating structure | [CRAFT Library](../CRAFT/CRAFT-Library.md) |
| Build a simulator for your own vessel | [How to Build Your Own Marine Craft Simulator](BUILD-YOUR-OWN-SIMULATOR.md) |
| Find guidance laws, controllers, or control allocation functions | [GNC Library](../GNC/GNC-Library.md) |
| Try attitude estimation or aided inertial navigation | [INS Library](../INS/INS-Library.md) |
| Work through a textbook example | [Examples in the MSS Quick Reference](../MSS%20Quick%20Reference.md#mss-examples-m-files) |

Run `mssHelp` to open the online [MSS Quick Reference](../MSS%20Quick%20Reference.md).
For the theory behind the models and algorithms, see the
[handbook lecture notes](https://wiley.fossen.biz/) and [MSS references](REFERENCES.md).
