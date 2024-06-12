function pathSignals = getPathSignals(path,s)
% psthSignals is compatibel with MATLAB and GNU Octave (www.octave.org).
% pathSignals = getPathSignals(path,s) generates the coefficients for 
% subpaths between given waypoints. Each subpath is parametrized by t in 
% [0,1), and smoothly connected at the waypoints. This corresponds to a 
% hybrid parametrized path, where Index i identifies the subpath, and the 
% parameter t in [0,1)identifies the location along the subpath.
%
% Input data:
%  path         - A Matlab data structure with the following fields:
%                .NumSubpaths: Number of subpaths.
%                .Order: Order of polynomials.
%                .WP: x- and y-coordinates of waypoints.
%                .LinSys: Liner set of equations Ax=b to solve the subpaths:
%                   .A:     Common A matrix for both x- and y-coefficients.
%                   .bx:    Cell structure that contains the b-vector for 
%                           each subpath for the x-coordinates.
%                   .by:    Cell structure that contains the b-vector for 
%                           each subpath for the y-coordinates.
%                .coeff: Coefficients for the subpaths:
%                   .a:     Cell structure that contains the a-coefficients 
%                           for the x-coordinates.
%                   .b:     Cell structure that contains the b-coefficients 
%                           for the y-coordinates.
%                   .a_der: Cell structure that contains the a-coefficients 
%                           for the derivatives 
%                           for the x-coordinates of the path.
%                   .b_der: Cell structure that contains the b-coefficients 
%                           for the derivatives for the y-coordinates of 
%                           the path.
%               - Scalar path variable in interval [0,N).
%
% Output data:
%  pathSignals  - A Matlab data structure with the following fields:
%                .pd:     Location (xd,yd) along the path corresponding to 
%                         the value s.
%                .pd_der: Value of path derivatives at (xd,yd) 
%                         corresponding to s. This contains a cell for each 
%                         derivative vector up to the polynomial order of 
%                         the path basis functions, i.e., 
%                         pathSignals.pd_der{1} contains 1st derivatives, 
%                         pathSignals.pd_der{2} the 2nd derivatives, and so on.
%                         Note that the path derivatives are guaranteed to 
%                         be continuous at the waypoints only for 
%                         derivatives up to r = (Order-1)/2.
%
%    Author:        Roger Skjetne
%    Date created:  2014-02-10  Roger Skjetne
%    Revised:      	

%% Initialization
pathSignals = [];
ord = path.Order;
if s < 0
    s = 0.0;
elseif s >= path.NumSubpaths
    s = path.NumSubpaths-1e3*eps;
end
idx = floor(s)+1;
t   = s-(idx-1);


%% Getting values
xd = 0;
yd = 0;
for k=0:ord
    xd = xd + path.coeff.a{idx}(k+1)*t^(ord-k);
    yd = yd + path.coeff.b{idx}(k+1)*t^(ord-k);
end
pathSignals.pd = [xd yd];
for k=1:ord
    xd = 0; yd = 0;
    a_vec = fliplr(path.coeff.a_der{idx}{k});
    b_vec = fliplr(path.coeff.b_der{idx}{k});
    for j=1:length(a_vec)
        xd = xd + a_vec(j)*t^(j-1);
    end
    for j=1:length(b_vec)
        yd = yd + b_vec(j)*t^(j-1);
    end
    pathSignals.pd_der{k} = [xd yd];
end
    
end