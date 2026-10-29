function drift = readWamitMeanDrift(filename, rho, g, ULEN)
% readWamitMeanDrift reads signed diagonal mean-drift coefficients from .8/.9.
% Supporting reader used by wamit2vessel.
% drift = readWamitMeanDrift(filename,rho,g,ULEN) returns the current MSS
% drift.w and drift.coeff{1:3} layout for physical DOFs 1, 2 and 6.
% Coefficients are loads per squared incident-wave amplitude: N/m^2 for
% surge/sway and Nm/m^2 for yaw. The reference remains the WAMIT body origin;
% this routine does not translate moments to CG. Axes are MSS FSD.
% Input must contain diagonal data on the 0:10:180 degree heading grid.
% The opposite half-plane is filled using port-starboard symmetry.
%
% Author: Thor I. Fossen
% Date:   2026-10-28

validateattributes(rho, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(g, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(ULEN, {'numeric'}, {'scalar','real','finite','positive'});
% Some MAT-files store density as an integer; scaling must use real arithmetic.
rho = double(rho); g = double(g); ULEN = double(ULEN);
data = load(filename, '-ascii');
if size(data,2) ~= 8 || isempty(data) || any(~isfinite(data(:))) || ...
        any(data(:,1) <= 0) || any(abs(data(:,2)-data(:,3)) > 1e-4)
    error('MSS:WamitMeanDrift', ...
        'Expected finite diagonal WAMIT .8/.9 data with positive periods.');
end
headings = unique(data(:,2));
if ~isequal(headings(:), (0:10:180)')
    error('MSS:WamitMeanDrift', 'WAMIT drift headings must be 0:10:180 degrees.');
end
[omega, order] = sort(2*pi./unique(data(:,1)));
periods = unique(data(:,1));
periods = periods(order);
drift.w = omega(:)';
drift.coeff = cell(1,3);
physicalDOFs = [1 2 6];
% Same axes/heading conversion as the first-order WAMIT importer.
signs = [1 1 -1];
scales = rho*g*[ULEN ULEN ULEN^2];
for j = 1:3
    selected = data(data(:,4) == physicalDOFs(j),:);
    [~, freqIndex] = ismember(selected(:,1), periods);
    [~, headingIndex] = ismember(selected(:,2), headings);
    linearIndex = sub2ind([numel(omega),19], freqIndex, headingIndex);
    if numel(linearIndex) ~= numel(omega)*19 || ...
            numel(unique(linearIndex)) ~= numel(omega)*19
        error('MSS:WamitMeanDrift', ...
            'Missing or duplicate diagonal drift data for DOF %d.', physicalDOFs(j));
    end
    values = zeros(numel(omega),36);
    values(linearIndex) = signs(j)*scales(j)*selected(:,7);
    parity = 1;
    if j > 1, parity = -1; end
    values(:,20:36) = parity*values(:,18:-1:2);
    drift.coeff{j} = values;
end
drift.reference = 'WAMIT body origin';
drift.axes = 'MSS FSD';
end
