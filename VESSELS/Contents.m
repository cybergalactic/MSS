% The MSS marine craft models are compatibel with MATLAB and the free 
% software GNU Octave (www.octave.org):
%
% clarke83        - Linear maneuvering model parametrized using (L,B,T) 
%                   found from linear regression of model tests (Clarke et 
%                   al. 1983).
% container       - Nonlinear maneuvering model of a high-speed container 
%                   ship, L = 175 m, including roll (Son and Nomoto 1982)
% DSRV            - Deep submergence rescue vehicle (DSRV), L = 5.0 m 
%                   (Healey 1992).
% frigate         - Nonlinear autopilot model for a frigate, L = 100 m.
% Lcontainer      - Linearized model of a high-speed container ship, 
%                   L = 175 m, including the roll mode (Son and Nomoto 1982).
% mariner         - Nonlinear maneuvering model for the Mariner class 
%                   vessel, L = 160 m. 
% navalvessel     - Nonlinear maneuvering model of a multipurpose naval 
%                   vessel, L = 51.5 m.
% npsauv          - Naval Postgraduate School autonomous underwater vehicle
%                   (AUV), L = 5.3 m. 
% osv             - Nonlinear model of an Offshore supply vessel (OSV), 
%                   L = 83.0 m.
% otter           - OTTER small autonomous USV, L = 2.0 m.
% remus100        - REMUS 100 autonomous underwater vehicle (AUV), L = 1.9 m
% rig             - Semisubmersible linear mass-damper-spring model, 
%                   L = 84.6 m.
% supply          - Linear DP model of a supply vessel, L = 76.2 m.
% tanker          - Nonlinear course unstable maneuvering model of a 
%                   tanker, L = 304 m. 
% zeefakkel       - Nonlinear autopilot model of a pleasure craft, L = 45 m.
%
% Script for time-series simulations:
%
% SIMclarke83     - Simulate clarke83.m under PD control.
% SIMdsrv         - Simulate DSRV.m with an autopilot for depth control 
%                   using successive-loop closure.
% SIMfrigate      - Simulate frigate.m using a PID heading autopilot.
% SIMmariner      - Simulate mariner.m under PD control.
% SIMotter        - Simulate otter.m under feedback control.
% SIMosv          - Simulate osv.m under nonlinear DP control with 
%                   constrained control allocation (dynamic optimization).
% SIMcontainer    - Simulate container.m and Lcontainer.m under PD control.
% SIMnavalvessel  - Simulate navalvessel.m under PD control.
% SIMnpsauv       - Simulate npsauv.m with MIMO PID autopilots for depth 
%                   and heading control, and 3-D straight-line path 
%                   following using ALOS
% SIMremus100     - Simulate remus100.m using autopilots for depth and 
%                   heading control, and adaptive line-of-sight (ALOS)
%                   guidance laws for 3-D path-following control.
% SIMrig          - Simulate the 6-DOF semisubmersible model under PID 
%                   control.
% SIMsupply       - Simulate the linear supply vessel model under DP control.
% SIMtanker       - Simulate the tanker model under DP control.