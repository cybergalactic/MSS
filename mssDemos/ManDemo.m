%       MANEUVERING TRIALS demonstrated by simulating m-file vessel models: 
%       Type: >> help gnc\vesselmodels  for more information
%
%       1)  Turning circle for the Mariner class vessel and a container ship 
%           (see mariner.m and container.m)
%       2)  Zig-zag manuvers for the Mariner class vessel and a container ship 
%           (see mariner.m and container.m)
%       3)  Pullout maneuver for the Mariner class vessel and the Esso Osaka tanker 
%           (see mariner.m and tanker.m)
%
%       0)  Quit
echo off

while 1
    demos= ['ExTurnCircle'
            'ExZigZag    '
            'ExPullout   '];
    clc
    help mandemo
    n = input('Select a demo number: ');
    if ((n <= 0) | (n > 3)) 
        break
    end
    demos = demos(n,:);
    eval(demos)
    clear
 end
 clear n demos
clc
