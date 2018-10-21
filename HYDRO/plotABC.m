function plotABC(vessel,mtrx,i,j,velno)
% PLOTABC (MSS Hydro)
%
% >> plotABC(vessel,mtrx,i,j,velno) plots added mass, damping, restoring 
%    matrix element i,j versus frequency ans speed
%
% >> plotABC(vessel,mtrx) plots all elements Aij,Bij,Cij versus
%    frequency and speed
%
% Inputs: 
%   vessel:  MSS vessel structure generated from a data file: *
%   (file name without extension) by using the following calls: 
%
%   ShipX (VERES), MARINTEK :  vessel = veres2vessel(*)
%   WAMIT                   :  vessel = wamit2vessel(*)
%
%   mtrx = 'A'    added mass
%          'B'    potential + viscous damping
%          'C'    restoring forces
%
%  i, j (optionally):    matrix element
%  velno (optionally):   speed number
%
% Author:    Thor I. Fossen
% Date:      2005-05-16 First version
% Revisions: 2005-09-23 Removed viscous terms when plotting
% ________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

if strcmp(mtrx,'A')
    H = vessel.A;
    figno = 10;
elseif strcmp(mtrx,'B')
    H  = vessel.B; 
    figno = 20;
elseif strcmp(mtrx,'C')
    H = vessel.C;
    figno = 30;
end

% PLOT DATA
freqs      = vessel.freqs;
velocities = vessel.velocities;

nvel  = length(velocities);
nfreq = length(freqs);

if nargin == 2
    
    % A, B, C plots
    
    for velno = 1:nvel
        
        % longitudinal plots
        k = 1;
        figure(figno)
        
        for i = 1:2:5
            
            for j = 1:2:5
                
                Hplot = reshape(H(i,j,:,velno),nfreq,1);
                splot = 330+k;
                subplot(splot)
                plot(freqs,Hplot(:,1),'b-o'),grid
                    
                if mtrx == 'A'
                    Hw = strcat(strcat(strcat(...
                        'A_{',num2str(i)),num2str(j)),'}');
                elseif mtrx == 'B'
                    Hw = strcat(strcat(strcat(...
                        'B_{',num2str(i)),num2str(j)),'}');
                elseif mtrx == 'C'
                    Hw = strcat(strcat(strcat(...
                        'C_{',num2str(i)),num2str(j)),'}');
                end
                U = strcat(strcat(...
                    ' (U=', num2str(velocities(velno))),' m/s)');
                
                xlabel('frequency (rad/s)')
                title(strcat(Hw,U));
                k = k +1;
                
            end
        end
        
        figno = figno + 1;
        
        % lateral plots
        k = 1;
        figure(figno)
        
        for i = 2:2:6
            
            for j = 2:2:6
                
                Hplot = reshape(H(i,j,:,velno),nfreq,1);
                splot = 330+k;
                subplot(splot)   
                plot(freqs,Hplot(:,1),'b-o'),grid
                
                if mtrx == 'A'
                    Hw = strcat(strcat(strcat(...
                        'A_{',num2str(i)),num2str(j)),'}');
                elseif mtrx == 'B'
                    Hw = strcat(strcat(strcat(...
                        'B_{',num2str(i)),num2str(j)),'}');
                elseif mtrx == 'C'
                    Hw = strcat(strcat(strcat(...
                        'C_{',num2str(i)),num2str(j)),'}');
                end
                U = strcat(strcat(...
                    ' (U=', num2str(velocities(velno))),' m/s)');
                
                xlabel('frequency (rad/s)')
                title(strcat(Hw,U));
                k = k +1;
                
            end
        end
        
        figno = figno + 1;
        
    end

elseif nargin == 5

    figure(3)
    subplot(111)
    Hplot = reshape(H(i,j,:,velno),nfreq,1);
    
    % A, B, C plots
    
    plot(freqs,Hplot(:,1),'b-o'),grid
    
    if mtrx == 'A'
        Hw = strcat(strcat(strcat(...
            'A_{',num2str(i)),num2str(j)),'}');
    elseif mtrx == 'B'
        Hw = strcat(strcat(strcat(...
            'B_{',num2str(i)),num2str(j)),'}');
    elseif mtrx == 'C'
        Hw = strcat(strcat(strcat(...
            'C_{',num2str(i)),num2str(j)),'}');
    end
    U = strcat(strcat(...
        ' (U=', num2str(velocities(velno))),' m/s)');
    
    xlabel('frequency (rad/s)')
    title(strcat(Hw,U));

else
    disp('Error: wrong number of input arguments')
    return
end
