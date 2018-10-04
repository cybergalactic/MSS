function [table, pathlength] = order5(drawmode, constraint, x_vector,y_vector)
% [table, pathlength] = ORDER5(drawmode, constraint, x_vector,y_vector)
% Path generation using 5th order polynominals. To be used with demowpt.mdl
%
% This version uses a fifth-order polynomial and makes the second
% derivatives in each waypoint start and end in zero. The calculations are
% done offline, and the coeficcients are returned in table-form. 
%
% Example call: [table, pathlength] = order5(1, 0, [0 200 400 700 1000],[0 600 500 400 1200])
%
% Author:       Jørgen Corneliussen
% Date:         20 Nov 2002
% Revisions: 


% Testing the provided waypoints
if length(x_vector) ~= length(y_vector) %There is not an equal amount of x- and y-coordinates.
    disp('The number of X-coordinates and Y-coordinates must be the same');
elseif length(x_vector) == 1 % Only one waypoint provided
    disp('You must provide more than one waypoint');
elseif length(x_vector) == 2 %A straight line, inserting another waypoint in the middle
    disp('Only two waypoints were provided. Adding another waypoint in the middle');
    x_vector = [x_vector(1) (x_vector(2) - x_vector(1))/2 x_vector(2)];
    y_vector = [y_vector(1) (y_vector(2) - y_vector(1))/2 y_vector(2)];
end

N = length(x_vector)-1; %The number of paths

dx = diff(x_vector); % Finding the desired derivatives in each waypoint
dy = diff(y_vector);

dx(N+1) = 0;%Setting the derivatives in the last waypoint to zero.
dy(N+1) = 0;

if constraint % If large differences in the derivatives occur, decrease them
    for i=1:N
        if (abs(dx(i+1)) >= abs(3*dx(i)) & dx(i) ~= 0), dx(i+1) = 3*dx(i); end
        if (abs(dy(i+1)) >= abs(3*dy(i)) & dy(i) ~= 0), dy(i+1) = 3*dy(i); end
    end
end

ddx = zeros(1,N+1); % Setting the second derivatives to zero in all waypoints
ddy = zeros(1,N+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% COMPUTATION OF THE PATH ( x(th),y(th) ) %%%%%%%%%%%%%%%%%
%%%%%  x(th) = a5*th^5 + a4*th^4 + a3*th^3 + a2*th^2 + a1*th + a0  %%%%%
%%%%%  y(th) = b5*th^5 + b4*th^4 + b3*th^3 + b2*th^2 + b1*th + b0  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:N,
    
    % Select two active way-points for curve fitting
    x0     = x_vector(i);
    x1     = x_vector(i+1);    
    y0     = y_vector(i);
    y1     = y_vector(i+1);
    th0    = 0;
    th1    = 1;
    
    % y = A*x,   x = A'*inv(A*A')*y
    % unknown parameters: x = [a5 a4 a3 a2 a1 a0 b5 b4 b3 b2 b1 b0]'
    
  y = [x0 x1 y0 y1 dx(i) dx(i+1) dy(i) dy(i+1) ddx(i) ddx(i+1) ddy(i) ddy(i+1)]';
  
  A = [ th0^5       th0^4       th0^3       th0^2       th0         1           0           0           0           0       0       0        % x(th0) = x0
        th1^5       th1^4       th1^3       th1^2       th1         1           0           0           0           0       0       0        % x(th1) = x1
        0           0           0           0           0           0           th0^5       th0^4       th0^3       th0^2   th0     1        % y(th0) = y0
        0           0           0           0           0           0           th1^5       th1^4       th1^3       th1^2   th1     1        % y(th1) = y1  
        5*th0^4     4*th0^3     3*th0^2     2*th0       1           0           0           0           0           0       0       0
        5*th1^4     4*th1^3     3*th1^2     2*th1       1           0           0           0           0           0       0       0
        0           0           0           0           0           0           5*th0^4     4*th0^3     3*th0^2     2*th0   1       0
        0           0           0           0           0           0           5*th1^4     4*th1^3     3*th1^2     2*th1   1       0
        20*th0^3    12*th0^2    6*th0       2           0           0           0           0           0           0       0       0
        20*th1^3    12*th1^2    6*th1       2           0           0           0           0           0           0       0       0
        0           0           0           0           0           0           20*th0^3    12*th0^2    6*th0       2       0       0
        0           0           0           0           0           0           20*th1^3    12*th1^2    6*th1       2       0       0];

   x = inv(A)*y;

   Acoeff(i,:) = x(1:6)';  %[a5 a4 a3 a2 a1 a0]
   Bcoeff(i,:) = x(7:12)'; %[b5 b4 b3 b2 b1 b0]
end

pathlength = 0;
for i= 1:N, %Calculating the length of the path
    th = linspace(0,1);
    xpath = polyval(Acoeff(i,:),th);
    ypath = polyval(Bcoeff(i,:),th);
    for j=1:length(xpath)-1
        pathlength = pathlength + sqrt( (xpath(j+1)-xpath(j))^2 + (ypath(j+1)-ypath(j))^2);
    end    
end

if drawmode == 1 %Plotting the path if drawmode is set to 1
    clf
    hold on
    plot(x_vector,y_vector,'ro', 'linewidth',2)
    for i= 1:N,
        th = linspace(0,1);
        plot(polyval(Acoeff(i,:),th),polyval(Bcoeff(i,:),th),'linewidth',2);
    end
    grid
    title('The desired path computed with a 5th order algorithm')
    xlabel('X-direction')
    ylabel('Y-direction')
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                CREATING THE DERIVATIVES OF THE PATH                          %%
%%         The first derivatives by means of theta are:                         %%
%%           x = 5*a5*theta^4 + 4*a4*theta^3 + 3*a3*theta^2 + 2*a2*theta + a1   %%
%%           y = 5*b5*theta^4 + 4*b4*theta^3 + 3*b3*theta^2 + 2*b2*theta + b1   %%
%%         The second derivatives by means of theta are:                        %%
%%                x = 20*a5*theta^3 + 12*a4^2 + 6*a3*theta + 2*a2               %%
%%                y = 20*b5*theta^3 + 12*b4^2 + 6*b3*theta + 2*b2               %%
%%         The third derivatives are:                                           %%
%%                x = 60*a5*theta^2 + 24*a4 + 6*a3                              %%
%%                y = 60*b5*theta^2 + 24*b4 + 6*b3                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:N
    x_theta_one(i,:) = [0 5*Acoeff(i,1) 4*Acoeff(i,2) 3*Acoeff(i,3) 2*Acoeff(i,4) Acoeff(i,5)];
    y_theta_one(i,:) = [0 5*Bcoeff(i,1) 4*Bcoeff(i,2) 3*Bcoeff(i,3) 2*Bcoeff(i,4) Bcoeff(i,5)];
    
    x_theta_two(i,:) = [0 0 20*Acoeff(i,1)  12*Acoeff(i,2) 6*Acoeff(i,3) 2*Acoeff(i,4)];
    y_theta_two(i,:) = [0 0 20*Bcoeff(i,1)  12*Bcoeff(i,2) 6*Bcoeff(i,3) 2*Bcoeff(i,4)];
    
    x_theta_three(i,:) = [0 0 0 60*Acoeff(i,1)  24*Acoeff(i,2) 6*Acoeff(i,3)];
    y_theta_three(i,:) = [0 0 0 60*Bcoeff(i,1)  24*Bcoeff(i,2) 6*Bcoeff(i,3)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CREATING THE LOOKUP TABLE CONTAINING ALL INFORMATION ABOUT THE PATH         %%
%%   This table is for use in demowpt.mdl                                       %%
%%   The table has N-1 columns, where N is the number of waypoints.             %%
%%                                                                              %%
%%  - Column 1-6 contains the polynomial for the x-direction                    %%
%%  - Column 7-12 contains the polynomial for the y-direction                   %%
%%  - Column 13-18 contains the polynomial for x differentiated once with       %%
%%      regards to theta                                                        %%
%%  - Column 19-24 contains the polynomial for y differentiated once with       %%
%%      regards to theta                                                        %%
%%  - Column 25-30 contains the polynomial for x differentiated twice           %%
%%      with regards to theta                                                   %%
%%  - Column 31-36 contains the polynomial for x differentiated twice           %%
%%      with regards to theta                                                   %%
%%  - Column 37-42 contains the polynomial for x differentiated three           %%
%%      times with regards to theta                                             %%
%%  - Column 43-48 contains the polynomial for x differentiated three           %%
%%      times with regards to theta                                             %%
%%                                                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table(:,1:6) = Acoeff; %Inserting the path polynomials into the table. 
table(:,7:12) = Bcoeff; %Inserting the path polynomials into the table.

table(:,13:18) = x_theta_one;  
table(:,19:24) = y_theta_one; 

table(:,25:30) = x_theta_two; 
table(:,31:36) = y_theta_two; 

table(:,37:42) = x_theta_three; 
table(:,43:48) = y_theta_three; 

table = table';