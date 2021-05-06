% MSS Vessel Models
%
% Standard ship, semisubmersible and underwater vehicle models.
%
%    Ships:
%    clarke83    - Ship maneuvering model parametrized using L, B and T
%    frigate     - Frigate, L = 100 m (nonlinear autopilot model)
%    mariner     - Mariner class vessel, L=160 m (nonlinear maneuvering model)
%    tanker      - Tanker, L=304 m (nonlinear course unstable maneuvering model)
%    container   - Container ship, L=175 m (nonlinear maneuvering model including the roll mode)
%    Lcontainer  - Container ship, L=175 m (LINEAR maneuvering model including the roll mode)
%    navalvessel - Multipurpose naval vessel, L = 51.5 m (nonlinear manneuvering model)
%    ROVzefakkel - Boat, L = 45 m (nonlinear autopilot model) 
%    supply      - Supply vessel, L = 76.2 m (linear DP model)
%
%    Underwater vehicles:
%    DSRV        - Deep submergence rescue vehicle (DSRV), L = 5.0 m
%    npsauv      - Naval Postgraduate School autnomous underwater vehicle (AUV), L = 5.3 m
%	 remus100    - Remus 100 autnomous underwater vehicle (AUV), L = 1.9 m

%    Semisubmersibles:
%    rig         - Semisub (mass-damper-spring) model, L = 84.6 m
%
%    USV models:
%    otter       - Small autonomous USV, L = 2.0 m
%
%    Time-series simulations:
%    SIMclarke83    - Simulate clarke83.m under PD control
%    SIMmariner     - Simulate mariner.m under PD control
%    SIMotter       - Smulate otter.m under feedback control
%    SIMcontainer   - Simulate container.m and Lcontainer.m under PD control
%    SIMnavalvessel - Simulate navalvessel.m under PD control
%    SIMrig         - Simulate the 6-DOF semisub model under PID control
