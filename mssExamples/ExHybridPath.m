% ExHybridPath Function that calculates the hybrid continuous path based on waypoints.
%
%    Author:        Roger Skjetne
%    Date created:  2014-02-10  Roger Skjetne.
%    Revised: 

%% Initialization
X_wp = [0 5 10 15 20 25 30 35 40];      % Waypoints
Y_wp = [0 9 10 11 20 10 20 10 20];
% X_wp = [10  5  0  5 10  5  0  5 10  5  0  5 10];
% Y_wp = [10 15 10  5 10 15 10  5 10 15 10  5 10];

s           = 1.0;  % Path parameter in [0, NumSubPaths]
lambda      = 0.25; % Curvature constant.
r           = 2;    % Differentiability order.
PlotHandle  = 1;    % Empty if no plotting

%% Run scripts
WP = [X_wp' Y_wp'];
path = hybridPath(WP,r,lambda,PlotHandle);
pathSignals = getPathSignals(Path,s);

%% Plotting
if PlotHandle
    figure(PlotHandle); hold on;
    plot(pathSignals.pd(1),pathSignals.pd(2),'sm','LineWidth',1.5);
    plot([pathSignals.pd(1) pathSignals.pd(1)+pathSignals.pd_der{1}(1)],...
    [pathSignals.pd(2) pathSignals.pd(2)+pathSignals.pd_der{1}(2)],'m','LineWidth',1.25);

    s = [0:.01:path.NumSubpaths];
    for k=1:Path.Order
        pd_der = zeros(length(s),2);
        for j=1:length(s)
            PS = getPathSignals(Path,s(j));
            pd_der(j,:) = PS.pd_der{k};
        end
    
        % Plotting derivatives
        figure(PlotHandle+k); clf; 
        plot(s,pd_der(:,1),s,pd_der(:,2),'LineWidth',1.25); 
        text_x = ['x_d^{s^',num2str(k),'}'];
        text_y = ['y_d^{s^',num2str(k),'}'];
        legend(text_x,text_y); grid on;
    end
end

%% Function 
function pathSignals = getPathSignals(path,s)
% pathSignals = getPathSignals(path,s)
%
% Function that generates the coefficients for subpaths between given waypoints.
% Each subpath is parametrized by t in [0,1), and smoothly connected at the waypoints. This 
% then correspond to a hybrid parametrized path, where Index i identifies the subpath, and 
% the parameter t in [0,1) identifies the location along the subpath.
%
% Input data:
%  path         - A Matlab data structure with the following fields:
%                .NumSubpaths: Number of subpaths.
%                .Order: Order of polynomials.
%                .WP: x- and y-coordinates of waypoints.
%                .LinSys: Liner set of equations Ax=b to solve the subpaths:
%                   .A:     Common A matrix for both x- and y-coefficients.
%                   .bx:    Cell structure that contains the b-vector for each subpath for the x-coordinates.
%                   .by:    Cell structure that contains the b-vector for each subpath for the y-coordinates.
%                .coeff: Coefficients for the subpaths:
%                   .a:     Cell structure that contains the a-coefficients for the x-coordinates.
%                   .b:     Cell structure that contains the b-coefficients for the y-coordinates.
%                   .a_der: Cell structure that contains the a-coefficients for the derivatives 
%                           for the x-coordinates of the path.
%                   .b_der: Cell structure that contains the b-coefficients for the derivatives 
%                           for the y-coordinates of the path.
%  s            - Scalar path variable in interval [0,N).
%
% Output data:
%  pathSignals  - A Matlab data structure with the following fields:
%                .pd: Location (xd,yd) along the path corresponding to the value s.
%                .pd_der: Value of path derivatives at (xd,yd) corresponding to s. 
%                         This contains a cell for each derivative vector up to the 
%                         polynomial order of the path basis functions, i.e., 
%                         pathSignals.pd_der{1} contains 1st derivatives, 
%                         pathSignals.pd_der{2} the 2nd derivatives, and so on.
%                 Note that the path derivatives are guaranteed to be continuous 
%                 at the waypoints only for derivatives up to r = (Order-1)/2.
%
%    Author:        Roger Skjetne
%    Date created:  2014-02-10  Roger Skjetne.
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


