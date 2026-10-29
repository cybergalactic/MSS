% Hydrodynamic processing of WAMIT semi-sub data
vessel = wamit2vessel('fpso',12,200,44);
vesselABC = vessel2ss(vessel,20,5,8,1);

% plot
plotABC(vessel,'A')
plotABC(vessel,'B')
plotTF(vessel,'motion','rads',1)
plotTF(vessel,'force','rads',1)
plotWD(vessel,'rads',1)

% display main data
display(vessel.main);

% plot WAMIT geometry file
plot_wamitgdf('fpso_low');


