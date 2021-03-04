%       ------- GUIDANCE, NAVIGATION AND CONTROL TOOLBOX Demonstrations---
%
%       1)  Euler angle and quaternion kinematics
%       2)  Ship maneuvering trials 
%       3)  User editable script for simulation of m-file vessel models
%       4)  Straight-line, directional and positional motion stability
%       5)  Wave spectrum demonstration
%
%       0)  Quit
echo off

%   Guidance, Navigation and Control Toolbox Demonstrations.
%
%   Thor I. Fossen 2001-08-14


while 1
    demos= ['KinDemo '
            'ManDemo '
            'SimDemo ' 
            'StabDemo'
            'WaveDemo'];
    clc
    help gncdemo
    n = input('Select a demo number: ');
    if ((n <= 0) | (n > 5)) 
        break
    end
    demos = demos(n,:);
    eval(demos)
    clear
 end
 clear n demos
clc
