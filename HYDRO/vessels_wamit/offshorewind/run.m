% run.m runs WAMIT and processes the WAMIT data file cylinder.out

% dos('c:\wamitv64\wamit.exe');    % REQUIRES WAMIT HARDWARE KEY

% Hydrodynamic processing of WAMIT semi-sub data
%vessel = wamit2vessel('cylinder',100,10,10);   % local version for cylindeer
load cylinder

% plot
plotABC(vessel,'A')
plotABC(vessel,'B')
%plotTF(vessel,'motion','rads',1)
plotTF(vessel,'force','rads',1)
plotWD(vessel,'rads',1)

% display main data
display(vessel.main);

% plot WAMIT geometry file (half cylinder)
plot_wamitgdf('cylinder_low');



