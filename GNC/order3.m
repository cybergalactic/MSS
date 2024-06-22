function [table, pathlength] = order3(x_vector,y_vector, drawmode, mode, xstart, xfinal, ystart, yfinal)
% [table, pathlength] = ORDER3(x_vector,y_vector, drawmode, mode, xstart, xfinal, ystart, yfinal)
%
% Path generation using cubic polynominals. To be used with demowpt.mdl
%
% This version uses the cubic spline method to create a third-order polynomial going through 
% all the waypoints. All calculations are made offline, and the coeficcients are returned in
% table-form. The variable is theta, and it is increased with one between each waypoint.
%
% Usage: 
% [table, pathlength] = order3([0 200 400 700 1000],[0 600 500 400 1200],1,1,0,0,0,0)
%
% Author:   Jørgen Corneliussen
% Date:     20 November 2002
% Revisions:

% Testing the provided waypoints
if length(x_vector) ~= length(y_vector) %There is not an equal amount of x- and y-coordinates.
    disp('The number of X-coordinates and Y-coordinates must be the same');
elseif isscalar(x_vector) % Only one waypoint provided
    disp('You must provide more than one waypoint');
elseif length(x_vector) == 2 %A straight line, inserting another waypoint in the middle
    disp('Only two waypoints were provided. Adding another waypoint in the middle');
    x_vector = [x_vector(1) (x_vector(2) - x_vector(1))/2 x_vector(2)];
    y_vector = [y_vector(1) (y_vector(2) - y_vector(1))/2 y_vector(2)];
end

N = length(x_vector)-1;     % The number of paths
th = 0:N;                   % The theta-values along the path

if mode == 1 % The velocity is set in the endpoints
    c1 = [3*th(1)^2 2*th(1) 1 0];
    c2 = [3*th(N+1)^2 2*th(N+1) 1 0];
else % The acceleration is set in the endpoints
    c1 = [6*th(1) 2 0  0];
    c2 = [6*th(N+1) 2 0  0]; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%        COMPUTATION OF THE PATH ( x(th),y(th) )     %%%%%%%%%%%
%%%%%%%%%        x(th) = a3*th^3 + a2*th^2 + a1*th + a0      %%%%%%%%%%%
%%%%%%%%%        y(th) = b3*th^3 + b2*th^2 + b1*th + b0      %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creation of the y-vector:

y_x = [xstart x_vector(1)];
y_y = [ystart y_vector(1)];
for i=2:N
    y_x = [y_x x_vector(i) x_vector(i) 0 0];
    y_y = [y_y y_vector(i) y_vector(i) 0 0];    
end
y_x = [y_x x_vector(N+1) xfinal]';
y_y = [y_y y_vector(N+1) yfinal]';


% y = A*x,   x = inv(A)*y
% Unknown parameters: x_x = [a3 a2 a1 a0]', y_y = [b3 b2 b1 b0]'.
% The A matrix for the calculation of the coeficcients. This matrix will be the same for both
% the x- and y-direction. See Ch. 5 in Fossen, T. I. (2002). Marine Control Systems: Guidance, 
% Navigation and Control of Ships, Rigs and Underwater Vehicles, Marine Cybernetics AS. 
% ISBN 82-92356-00-2, www.marinecybernetics.com

A = zeros(4*N,4*N); % Defining the A matrix as a matrix containing only zeros
A(1,1:4) = c1;      % Inserting the startvalue on either the velocity or the acceleration
A(2,1:4) = [th(1)^3 th(1)^2 th(1) 1];   % The position in the starting point
A(4*N,4*N-3:4*N) = c2;                  % The position in the endpoint
A(4*N-1,4*N-3:4*N) = [th(N+1)^3 th(N+1)^2 th(N+1) 1];  % Inserting the final value on either the velocity or the acceleration

for i=2:N % Filling in the other constraints.
    row = 4*(i-2)+3; % Selecting the right row in the A-matrix
    col = 4*(i-1)+1; % Selecting the right column in the A-matrix
    
    A(row, col-4:col-1)     = [ th(i)^3     th(i)^2     th(i)   1];
    
    A(row+1, col:col+3)     = [ th(i)^3     th(i)^2     th(i)   1];
    
    A(row+2, col-4:col-1)   =-[ 3*th(i)^2   2*th(i)     1       0];
    A(row+2, col:col+3)     = [ 3*th(i)^2   2*th(i)     1       0];    
    
    A(row+3, col-4:col-1)   =-[ 6*th(i)     2           0       0];
    A(row+3, col:col+3)     = [ 6*th(i)     2           0       0];    
end

x_x = inv(A)*y_x; % The actual computation of the coefficients for the graphs in x-direction and y-direction.
x_y = inv(A)*y_y; 

for i=1:N  % Selecting the coefficcients and finding the derivatives of the path
    Acoeff(i,:) = x_x(4*(i-1)+1 : 4*(i-1)+4)';  
    Bcoeff(i,:) = x_y(4*(i-1)+1 : 4*(i-1)+4)';

    x_theta_one(i,:) = [0 3*Acoeff(i,1) 2*Acoeff(i,2) Acoeff(i,3)];
    y_theta_one(i,:) = [0 3*Bcoeff(i,1) 2*Bcoeff(i,2) Bcoeff(i,3)];
    
    x_theta_two(i,:) = [0 0 6*Acoeff(i,1) 2*Acoeff(i,2)];
    y_theta_two(i,:) = [0 0 6*Bcoeff(i,1) 2*Bcoeff(i,2)];
    
    x_theta_three(i,:) = [0 0 0 6*Acoeff(i,1)];
    y_theta_three(i,:) = [0 0 0 6*Bcoeff(i,1)];
end

pathlength = 0;
for i= 1:N, %Calculating the length of the path
    th = linspace(i-1,i);
    xpath = polyval(Acoeff(i,:),th);
    ypath = polyval(Bcoeff(i,:),th);
    for j=1:length(xpath)-1
        pathlength = pathlength + sqrt( (xpath(j+1)-xpath(j))^2 + (ypath(j+1)-ypath(j))^2);
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CREATING THE LOOKUP TABLE CONTAINING ALL INFORMATION ABOUT THE PATH         %%
%%   This table is for use in demowpt.mdl                                       %%
%%   The table has N-1 columns, where N is the number of waypoints.             %%
%%                                                                              %%
%%  - Column 1-4 contains the polynomial for the x-direction                    %%
%%  - Column 5-8 contains the polynomial for the y-direction                    %%
%%  - Column 9-12 contains the polynomial for x differentiated once with        %%
%%      regards to theta                                                        %%
%%  - Column 13-16 contains the polynomial for y differentiated once with       %%
%%      regards to theta                                                        %%
%%  - Column 17-20 contains the polynomial for x differentiated twice           %%
%%      with regards to theta                                                   %%
%%  - Column 21-24 contains the polynomial for x differentiated twice           %%
%%      with regards to theta                                                   %%
%%  - Column 25-28 contains the polynomial for x differentiated three           %%
%%      times with regards to theta                                             %%
%%  - Column 29-32 contains the polynomial for x differentiated three           %%
%%      times with regards to theta                                             %%
%%                                                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table(:,1:4) = Acoeff; %Inserting the path polynomials into the table. 
table(:,5:8) = Bcoeff; %Inserting the path polynomials into the table. 

table(:,9:12) = x_theta_one;
table(:,13:16) = y_theta_one; 

table(:,17:20) = x_theta_two; 
table(:,21:24) = y_theta_two; 

table(:,25:28) = x_theta_three; 
table(:,29:32) = y_theta_three; 

table = table';

%% Drawing the path and the waypoints if drawmode is set to 1

if drawmode
    figure(1)
    clf
    plot(x_vector,y_vector,'ro', 'linewidth',2)
    hold on
    for i= 1:N,
        th = linspace(i-1,i);  
        plot(polyval(Acoeff(i,:),th),polyval(Bcoeff(i,:),th), 'linewidth', 2);
    end
    grid
    title('The desired path computed using cubic splines')
    xlabel('X-direction')
    ylabel('Y-direction')
    hold off
end