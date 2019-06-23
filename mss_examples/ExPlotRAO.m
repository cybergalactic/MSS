% ExPlotRAOs  Script for plotting motion and force RAOs
%
% Author:    Thor I. Fossen
% Date:      23 June 2019
% Revisions: 

load supply;

for DOF = 1:6
    figure(DOF);
    
    % Plot motion RAO
    w     = vessel.motionRAO.w;
    amp   = vessel.motionRAO.amp;
    phase = vessel.motionRAO.phase;
    subplot(411); plotRAOamp(w,amp,DOF);  
    subplot(412); plotRAOphs(w,phase,DOF);
    
    % Plot force RAO
    w     = vessel.forceRAO.w;
    amp   = vessel.forceRAO.amp;
    phase = vessel.forceRAO.phase;
    subplot(413); plotRAOamp(w,amp,DOF);
    subplot(414); plotRAOphs(w,phase,DOF);
end


function plotRAOamp(w,amp,DOF)
    velno = 1;   
    arg = amp{DOF}(:,:,velno);
    plot(w,arg(:,1),'-*k',w, arg(:,4),'-ko',w, arg(:,7),'-k<',w, arg(:,10),'-kx')
    xlabel('wave encounter frequency (rad/s)')
    legend('0 deg','30 deg','60 deg','90 deg')
    grid
end

function plotRAOphs(w,phase,DOF)
    velno = 1;    
    phs = (180/pi)*unwrap(phase{DOF}(:,:,velno));  
    plot(w,phs(:,1),'-*k',w, phs(:,4),'-ko',w, phs(:,7),'-k<',w, phs(:,10),'-kx')
    legend('0 deg','30 deg','60 deg','90 deg')
    xlabel('wave encounter frequency (rad/s)')
    grid
end


