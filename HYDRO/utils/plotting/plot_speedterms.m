function plot_speedterms(vessel,Anew,Bnew,Bv)
% >> plot_speedterms(vessel,Anew,Bnew,Bv) plots added mass and damping 
%    versus frequency and speed
%
% Author:    Thor I. Fossen
% Date:      2005-04-18
% Revisions: 2008-02-15  R1.0
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2008 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% E-mail: contact@marinecontrol.org
% URL:    <http://www.marinecontrol.org>

close all
warning off

% Raw A and B data since plot_speedterms.m are called before ABCtransform.m
% Final vessel.A and vessel.B are equal to Anew and Bnew
freqs      = vessel.freqs;
velocities = vessel.velocities;
A          = vessel.A;   
B          = vessel.B;  

nvel  = length(velocities);
nfreq = length(freqs);
 
if nvel == 1,
    U1 = strcat(strcat('U=', num2str(velocities(1))),' m/s');  
elseif nvel == 2,
    U1 = strcat(strcat('U=', num2str(velocities(1))),' m/s');    
    U2 = strcat(strcat('U=', num2str(velocities(2))),' m/s'); 
    legend_txt = {U1 U2};       
elseif nvel == 3,
    U1 = strcat(strcat('U=', num2str(velocities(1))),' m/s');    
    U2 = strcat(strcat('U=', num2str(velocities(2))),' m/s');
    U3 = strcat(strcat('U=', num2str(velocities(3))),' m/s');
    legend_txt = {U1 U2 U3};           
elseif nvel == 4,
    U1 = strcat(strcat('U=', num2str(velocities(1))),' m/s');    
    U2 = strcat(strcat('U=', num2str(velocities(2))),' m/s');
    U3 = strcat(strcat('U=', num2str(velocities(3))),' m/s');    
    U4 = strcat(strcat('U=', num2str(velocities(4))),' m/s');  
    legend_txt = {U1 U2 U3 U4};         
elseif nvel == 5,
    U1 = strcat(strcat('U=', num2str(velocities(1))),' m/s');    
    U2 = strcat(strcat('U=', num2str(velocities(2))),' m/s');
    U3 = strcat(strcat('U=', num2str(velocities(3))),' m/s');    
    U4 = strcat(strcat('U=', num2str(velocities(4))),' m/s');  
    U5 = strcat(strcat('U=', num2str(velocities(5))),' m/s');
    legend_txt = {U1 U2 U3 U4 U5};         
else
    disp('Error: the maximum number of velocities must be 5')
    return
end
   
pcolor = {'g-' 'b-' 'k-' 'r-' 'y-' 'r*'};
    
% zero speed terms
vessel_new.A     = Anew;
vessel_new.B     = Bnew;
vessel_new.Bv    = Bv;
vessel_new.freqs = freqs;
vessel_new.velocities = velocities;

plotABC(vessel_new,'A')
plotABC(vessel_new,'B')

if nvel > 1  
    
    figno = 1;
    % ---------------------------------------------------------------------
    % PLOT Anew elements
    % ---------------------------------------------------------------------
    H = Anew;

    % longitudinal plots
    k = 1;
    figure(figno)

    for i = 1:2:5,

        for j = 1:2:5,

            splot = 330+k;
            subplot(splot)
            hold on

            for velno = 1:nvel,

                Hplot = reshape(H(i,j,:,velno),nfreq,1);
                plot(freqs,Hplot(:,1),char(pcolor(velno)))

                if velno == nvel,
                    Hw = strcat(strcat(strcat(...
                        'A_{',num2str(i)),num2str(j)),'}');
                    legend(legend_txt);
                    grid
                end

            end

            hold off
            xlabel('frequency (rad/s)')
            title(Hw);

            k = k +1;

        end
    end

    figno = figno + 1;

    % lateral plots
    k = 1;
    figure(figno)

    for i = 2:2:6,

        for j = 2:2:6,

            splot = 330+k;
            subplot(splot)
            hold on

            for velno = 1:nvel,

                Hplot = reshape(H(i,j,:,velno),nfreq,1);
                plot(freqs,Hplot(:,1),char(pcolor(velno)))

                if velno == nvel,
                    Hw = strcat(strcat(strcat(...
                        'A_{',num2str(i)),num2str(j)),'}');
                    legend(legend_txt);
                    grid
                end

            end

            hold off
            xlabel('frequency (rad/s)')
            title(Hw);

            k = k +1;

        end

    end

    figno = figno + 1;

    
% ---------------------------------------------------------------------
% PLOT Bnew elements
% ---------------------------------------------------------------------
    H = Bnew;

    % longitudinal plots
    k = 1;
    figure(figno)

    for i = 1:2:5,

        for j = 1:2:5,

            splot = 330+k;
            subplot(splot)
            hold on

            for velno = 1:nvel,

                Hplot = reshape(H(i,j,:,velno),nfreq,1);
                plot(freqs,Hplot(:,1),char(pcolor(velno)))

                if velno == nvel,
                    Hw = strcat(strcat(strcat(...
                        'B_{',num2str(i)),num2str(j)),'}');
                    legend(legend_txt);
                    grid
                end

            end

            hold off
            xlabel('frequency (rad/s)')
            title(Hw);

            k = k +1;

        end
    end

    figno = figno + 1;

    % lateral plots
    k = 1;
    figure(figno)

    for i = 2:2:6,
        for j = 2:2:6,

            splot = 330+k;
            subplot(splot)
            hold on

            for velno = 1:nvel,
                
                if i==4 & j==4
                    B44 = reshape(Bnew(4,4,:,velno),nfreq,1);
                    if isfield(vessel,'roll') % only veres data
                        B44 = B44 + vessel.roll.Bv44(:,velno);
                    end
                    plot(freqs,B44,char(pcolor(velno)))
                else
                    Hplot = reshape(H(i,j,:,velno),nfreq,1);
                    plot(freqs,Hplot(:,1),char(pcolor(velno)))
                end
                
                if velno == nvel,
                    Hw = strcat(strcat(strcat(...
                        'B_{',num2str(i)),num2str(j)),'}');
                    legend(legend_txt);
                    grid
                end

            end

            hold off
            xlabel('frequency (rad/s)')
            
            if i==4 & j==4
                title('B_{44}+B_{44}^*')
            else
                title(Hw);
            end

            k = k +1;

        end
    end

end

warning on
