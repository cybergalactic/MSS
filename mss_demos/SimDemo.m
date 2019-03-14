%       -------  MSS TOOLBOX Demonstrations---
%
%       User editable script for simulation of the m-file vessel models: 
%       Type: >> help gnc\vesselmodels  for more information
%
%       1)  Simulation of the Mariner class vessel under PD-control 
%           (see mariner.m)
%       2)  Simulation of the linear and nonlinear container ship models under PD-control 
%           (see container.m and Lcontainer.m)
%
%       0)  Quit
echo off

%   Guidance, Navigation and Control Toolbox Demonstrations.
%
%   Thor I. Fossen 2001-07-25

while 1
    demos= ['SimDemo1'
            'SimDemo2'];
    clc
    help simdemo
    n = input('Select a demo number: ');
    if ((n <= 0) | (n > 2)) 
        break
    end
    demos = demos(n,:);
    eval(demos)
    clear
 end
 clear n demos
clc
