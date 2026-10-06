function displayHydroData(vessel)

fprintf('%s\n','------------------------------------------------------------------');
fprintf('%s\n','VESSEL MAIN CHARACTERISTICS');
fprintf('%s\n','------------------------------------------------------------------');
fprintf('%-50s %6.2f m \n', 'Length (L):', vessel.main.Lpp);
fprintf('%-50s %6.2f m \n', 'Beam (B):', vessel.main.B);
fprintf('%-50s %6.2f m \n', 'Draft (T):', vessel.main.T);
fprintf('%-50s %6.2f kg \n', 'Mass (m):', vessel.main.m);
fprintf('%-50s %6.2f kg/m^3 \n', 'Density of water (rho):', vessel.main.rho);
fprintf('%-50s %6.2f m^3 \n', 'Volume displacement (nabla):', vessel.main.nabla);
fprintf('%-50s [%2.1f %2.1f %2.1f] m \n', 'Center of gravity (r_bG):',...
    vessel.main.CG(1), vessel.main.CG(2), vessel.main.CG(3));
fprintf('%-50s [%2.1f %2.1f %2.1f] m \n', 'Center of buoyancy (r_bB):',...
    vessel.main.CB(1), vessel.main.CB(2), vessel.main.CB(3));
if isfield(vessel.main, 'CF')
    fprintf('%-50s [%2.1f %2.1f %2.1f] m \n', ...
        'Center of flotation (r_bF):', vessel.main.CF(1), ...
        vessel.main.CF(2), vessel.main.CF(3));
end
fprintf('%-50s %4.2f m \n', 'Transverse metacentric height (GM_T):', vessel.main.GM_T);
fprintf('%-50s %4.2f m \n', 'Longitudinal metacentric height (GM_L):', vessel.main.GM_L);

matrices = {'System inertia matrix: M = MRB + MA', vessel.M;...
    'Linear damping matrix: D', vessel.D; 'Restoring matrix: G', vessel.G};

for k = 1:size(matrices, 1)
    fprintf('%s\n','-----------------------------------------------------------------');
    fprintf('%-40s\n', matrices{k, 1});
    for i = 1:size(matrices{k, 2}, 1)
        for j = 1:size(matrices{k, 2}, 2)
            if matrices{k, 2}(i,j) == 0
                fprintf('         0 ');
            else
                fprintf('%10.2e ', matrices{k, 2}(i,j));
            end
        end
        fprintf('\n');
    end
end
