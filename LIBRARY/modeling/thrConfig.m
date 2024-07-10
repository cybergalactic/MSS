function T_thr = thrConfig(alpha,l_x,l_y)
% T_thr = thrConfig(alpha,lx,ly) computes the thruster configuration matrix
% given by tau = T_thr(alpha) * K_thr * u following the definitions of 
% Fossen (2021, Chapter 11)
%
% Outputs:
%   T_thr: Thruster configuration matrix, dimension 3 x length(alpha)
%
% Inputs:
%   alpha: Cell array with arguments, example ['T', 0.7, 'M'] defines the
%          first thruster 'T' as a tunnel thruster, the second
%          as an azimuth thruster with orientation 0.7 radians, and the 
%          third 'M' as a main propleller
%   l_x:   Vector of thruster/propeller x-coordinates w.r.t the CO
%   l_y:   Vector of thruster/propeller y-coordinates w.r.t the CO
%
% Example: Ship with two main propellers located at [-20 5] and [-20 -5],
% one azimuth thruster located at [10 0], and one bow thruster located at
% [15 0] with azimuth angle 0.1 rad. The thruster configuration matrix is
% computed by
%
% T_thr = thrConfig( {'M','M',0.1,'T'}, [-20, -20, 10, 15], [5, -5, 0, 0] )
% T_thr =
%    1.0000    1.0000    0.9950         0
%         0         0    0.0998    1.0000
%   -5.0000    5.0000    0.9983   15.0000
%
% Author:    Thor I. Fossen
% Date:      2024-03-17
% Revisions: 

T_thr = zeros(3,length(alpha));    % Initialize the output matrix
    
% Loop through the alpha vector and concatenate the column vectors
for i = 1:length(alpha)

    if isequal(alpha{i}, 'M')           % Main propeller

        T_thr(:,i) = [ 1 0 -l_y(i) ]';

    elseif isequal(alpha{i}, 'T')       % Tunnel thruster

        T_thr(:,i) = [ 0 1 l_x(i) ]';

    else                                % Azimuth thruster

        az = alpha{i};
        T_thr(:,i) = [ cos(az) sin(az)  l_x(i)*sin(az)-l_y(i)*cos(az) ]';

    end

end
