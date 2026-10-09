% CRAFT - Marine craft models and time-domain simulation scripts
%
% CRAFT provides established craft models and a unified workflow for
% user-defined 6-DOF vessels generated from ShipX, WAMIT, or Capytaine
% seakeeping data. Full documentation is available in CRAFT/CRAFT-Library.md.
%
% Hydrodynamic vessel workflow:
%   SIMhydroVessel  - Editable simulation template for a user-defined USV,
%                     AUV, ship, or floating structure. Includes irregular
%                     waves, force RAOs, currents, autopilot, and DP control.
%   hydroVessel     - Nonlinear 12-state, 6-DOF equations of motion using
%                     the common MSS vessel structure.
%
% AUV simulation scripts and models:
%   SIMdsrv         - DSRV depth-control simulation.
%   SIMnpsauv       - NPS AUV depth, heading, and 3-D path-following simulation.
%   SIMremus100     - REMUS 100 depth, heading, and 3-D path-following simulation.
%   DSRV            - Deep Submergence Rescue Vehicle model.
%   npsauv          - Naval Postgraduate School AUV model.
%   remus100        - REMUS 100 AUV model.
%
% USV simulation scripts and models:
%   SIMotter        - OTTER USV guidance and control simulation.
%   otter           - OTTER USV model.
%
% Ship and floating-structure simulation scripts:
%   SIMclarke83     - Generic ship maneuvering and heading-control simulation.
%   SIMcontainer    - Container and linear-container ship simulation.
%   SIMfrigate      - Frigate heading-autopilot simulation.
%   SIMmariner      - Mariner heading and path-following simulation.
%   SIMnavalvessel  - Multipurpose naval-vessel simulation.
%   SIMosv          - Offshore supply vessel DP and control-allocation simulation.
%   SIMsemisub      - Semisubmersible 6-DOF control simulation.
%   SIMsupply       - Linear supply-vessel DP simulation.
%   SIMtanker       - Course-unstable tanker heading-control simulation.
%   SIMzeefakkel    - Zeefakkel heading-autopilot simulation.
%
% Ship and floating-structure models:
%   clarke83        - Linear ship maneuvering model parameterized by L, B, and T.
%   container       - Nonlinear container-ship model with roll dynamics.
%   Lcontainer      - Linearized container-ship model with roll dynamics.
%   frigate         - Nonlinear frigate heading model.
%   mariner         - Nonlinear Mariner-class ship model.
%   navalvessel     - Nonlinear multipurpose naval-vessel model.
%   osv             - Nonlinear offshore supply vessel model.
%   rig             - Linear semisubmersible mass-spring-damper model.
%   supply          - Linear supply-vessel DP model.
%   tanker          - Nonlinear course-unstable tanker model.
%   zeefakkel       - Nonlinear Zeefakkel recreational-craft model.
%
% See also MSSHELP.
