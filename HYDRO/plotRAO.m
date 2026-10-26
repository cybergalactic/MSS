function plotRAO(vessel)
% Function for plotting motion and force RAOs.
%
%    load 'my_vessel' vessel
%    exPlotRAO(vessel)
%
% Author:    Thor I. Fossen
% Date:      2026-10-26
% Revisions: 


pre  = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
post = {' (m)',' (m)',' (m)',' (deg)',' (deg)',' (deg)'};
txtMa = strcat(pre,' Motion RAO - Amplitude',post); 
txtMp = strcat(pre,' Motion RAO - Phase (deg)'); 
post = {' (N)',' (N)',' (N)',' (Nm)',' (Nm)',' (Nm)'};
txtFa = strcat(pre,' Force RAO - Amplitude',post); 
txtFp = strcat(pre,' Force RAO - Phase (deg)'); 

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

    set(findall(gcf,'type','line'),'linewidth',1.5);
    set(findall(gcf,'type','text'),'FontSize',12);
    set(findall(gcf,'type','legend'),'FontSize',10);
end

end

%% FUNCTIONS
function plotRAOamp(w,amp,DOF)
    velno = 1;   
    arg = amp{DOF}(:,:,velno);
    S = 1; if (DOF > 3), S = 180/pi; end  % scale DOF 4,5,6 to degrees
    plot(w,S*arg(:,1),'-*',w, S*arg(:,4),'-o',w, S*arg(:,7),'-<',w, S*arg(:,10),'-x')
    xlabel('Frequency (rad/s)')
    legend('0 deg','30 deg','60 deg','90 deg')
    grid
end

function plotRAOphs(w,phase,DOF)
    velno = 1;    
    txt = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
    phs = (180/pi)*unwrap(phase{DOF}(:,:,velno));  
    plot(w,phs(:,1),'-*',w, phs(:,4),'-o',w, phs(:,7),'-<',w, phs(:,10),'-x')
    legend('0 deg','30 deg','60 deg','90 deg')
    xlabel('Frequency (rad/s)')
    grid
end


