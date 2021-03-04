% ExPlotRAO  Script for plotting motion and force RAOs
%
% Author:    Thor I. Fossen
% Date:      23 June 2019
% Revisions: 

load supply;

pre  = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
post = {' (m)',' (m)',' (m)',' (deg)',' (deg)',' (deg)'};
txtMa = strcat(pre,' motion RAO - amplitude',post); 
txtMp = strcat(pre,' motion RAO - phase (deg)'); 
post = {' (N)',' (N)',' (N)',' (Nm)',' (Nm)',' (Nm)'};
txtFa = strcat(pre,' force RAO - amplitude',post); 
txtFp = strcat(pre,' force RAO - phase (deg)'); 

for DOF = 1:6
    figure(DOF);
    
    % Plot motion RAOs
    w     = vessel.motionRAO.w;
    amp   = vessel.motionRAO.amp;
    phase = vessel.motionRAO.phase;
    subplot(411); plotRAOamp(w,amp,DOF);   title(txtMa(DOF));
    subplot(412); plotRAOphs(w,phase,DOF); title(txtMp(DOF));    
    
    % Plot force RAOs
    w     = vessel.forceRAO.w;
    amp   = vessel.forceRAO.amp;
    phase = vessel.forceRAO.phase;    
    subplot(413); plotRAOamp(w,amp,DOF);   title(txtFa(DOF));    
    subplot(414); plotRAOphs(w,phase,DOF); title(txtFp(DOF));     
end


function plotRAOamp(w,amp,DOF)
    velno = 1;   
    arg = amp{DOF}(:,:,velno);
    S = 1; if (DOF > 3), S = 180/pi; end  % scale DOF 4,5,6 to degrees
    plot(w,S*arg(:,1),'-*k',w, S*arg(:,4),'-ko',w, S*arg(:,7),'-k<',w, S*arg(:,10),'-kx')
    xlabel('frequency (rad/s)')
    legend('0 deg','30 deg','60 deg','90 deg')
    grid
end

function plotRAOphs(w,phase,DOF)
    velno = 1;    
    txt = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
    phs = (180/pi)*unwrap(phase{DOF}(:,:,velno));  
    plot(w,phs(:,1),'-*k',w, phs(:,4),'-ko',w, phs(:,7),'-k<',w, phs(:,10),'-kx')
    legend('0 deg','30 deg','60 deg','90 deg')
    xlabel('frequency (rad/s)')
    grid
end


