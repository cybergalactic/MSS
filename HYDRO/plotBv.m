function plotBv(vessel)
% plotBv plots the zero-speed potential, viscous and total damping 
%    B_total(w) = B(w) + Bv(w) as a function of the frequency: 
%
%    plotBv(vessel)
%
% Input: 
%    vessel:  MSS vessel structure 
%
% Author:    Thor I. Fossen
% Date:      2020-03-08 First version
% Revisions: 

w       = vessel.freqs;
Nfreq   = length(w);
velno   = 1;

B       = vessel.B;
Bv      = vessel.Bv;

figno = 50;

% Longitudinal plots
k = 1;
figure(figno)
        
for i = 1:2:5
    for j = 1:2:5
        Bplot = reshape(B(i,j,:,velno),Nfreq,1);
        Bvplot = reshape(Bv(i,j,:),Nfreq,1);
        splot = 330+k;
        subplot(splot)
        plot(w,Bplot+Bvplot,'k-','linewidth',2)
        hold on
        plot(w,Bplot(:,1),'b-o')
        plot(w,Bvplot,'r-','linewidth',2)   
        hold off
        grid       
        Hw = strcat(strcat(strcat('B_{',num2str(i)),num2str(j)),'}');          
        xlabel('frequency (rad/s)')
        title(Hw);
        k = k +1;
        
    end
    
end

% Lateral plots
figno = figno + 1;
k = 1;
figure(figno)

for i = 2:2:6  
    for j = 2:2:6       
        Bplot = reshape(B(i,j,:,velno),Nfreq,1);
        Bvplot = reshape(Bv(i,j,:),Nfreq,1);
        splot = 330+k;
        subplot(splot)
        plot(w,Bplot+Bvplot,'k-','linewidth',2)
        hold on
        plot(w,Bplot(:,1),'b-o')
        plot(w,Bvplot,'r-','linewidth',2)   
        hold off
        grid      
        Hw = strcat(strcat(strcat('B_{',num2str(i)),num2str(j)),'}');
        xlabel('frequency (rad/s)')
        title(Hw);
        k = k +1;        
    end
end
