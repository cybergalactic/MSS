% The MSS craft models are compatible with MATLAB and the free
% software GNU Octave (www.octave.org):
%
% clarke83        - Linear ship maneuvering model parameterized by the main
%                   dimensions (L, B, T) using regression formulas based on
%                   model-test data (Clarke et al., 1983).
% container       - Nonlinear maneuvering model of a high-speed container
%                   ship, L = 175 m, including roll dynamics
%                   (Son and Nomoto, 1982).
% DSRV            - Deep Submergence Rescue Vehicle (DSRV), L = 5.0 m
%                   (Healey, 1992).
% frigate         - Nonlinear heading-autopilot model of a frigate,
%                   L = 100 m.
% Lcontainer      - Linearized model of a high-speed container ship,
%                   L = 175 m, including roll dynamics
%                   (Son and Nomoto, 1982).
% mariner         - Nonlinear maneuvering model of the Mariner-class ship,
%                   L = 160 m.
% navalvessel     - Nonlinear maneuvering model of a multipurpose naval
%                   vessel, L = 51.5 m.
% npsauv          - Naval Postgraduate School autonomous underwater vehicle
%                   (AUV), L = 5.3 m.
% osv             - Nonlinear model of an offshore supply vessel (OSV),
%                   L = 83.0 m.
% otter           - OTTER autonomous surface vehicle (USV), L = 2.0 m.
% remus100        - REMUS 100 autonomous underwater vehicle (AUV),
%                   L = 1.9 m.
% rig             - Linear mass-spring-damper model of a semisubmersible
%                   offshore platform, L = 84.6 m.
% supply          - Linear dynamic-positioning (DP) model of a supply
%                   vessel, L = 76.2 m.
% tanker          - Nonlinear course-unstable tanker model,
%                   L = 304 m.
% zeefakkel       - Nonlinear autopilot model of the Zeefakkel pleasure
%                   craft, L = 13.7 m (45 ft).
%
% Scripts for time-domain simulations:
%
% SIMclarke83     - Simulate clarke83.m under PD heading control.
% SIMcontainer    - Simulate container.m and Lcontainer.m under PD control.
% SIMdsrv         - Simulate DSRV.m using successive-loop-closure depth
%                   control.
% SIMfrigate      - Simulate frigate.m using a PID heading autopilot.
% SIMmariner      - Simulate mariner.m under PD heading control.
% SIMnavalvessel  - Simulate navalvessel.m under PD heading control.
% SIMnpsauv       - Simulate npsauv.m using MIMO PID autopilots for depth
%                   and heading control, and ALOS guidance for 3-D
%                   straight-line path following.
% SIMosv          - Simulate osv.m under nonlinear DP control with
%                   constrained control allocation.
% SIMotter        - Simulate otter.m under feedback control.
% SIMremus100     - Simulate remus100.m using depth and heading autopilots,
%                   and adaptive line-of-sight (ALOS) guidance for
%                   3-D path following.
% SIMrig          - Simulate the 6-DOF semisubmersible model under
%                   PID control.
% SIMsupply       - Simulate the linear supply-vessel model under
%                   DP control.
% SIMtanker       - Simulate the tanker model under heading control.