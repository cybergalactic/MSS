function plot_speedterms(vessel,Anew,Bnew,Bv)
% >> plot_speedterms(vessel,Anew,Bnew,Bv) plots added mass and damping 
%    versus frequency and speed
%
% Author:    Thor I. Fossen
% Date:      2005-04-18
% Revisions: 
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

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
