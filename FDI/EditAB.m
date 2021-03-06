function [Amod,Bmod,Wmod]=EditAB(A,B,W)
% Function to prepare the data for identification: 
% Select frequency range and eliminate wildpoints.
%
% Use: [Amod,Bmod,Wmod] = EditAB(A,B,W)
%
% A - vector of frequency dependant added mass A(w)
% B - vector of frequency dependant potential damping B(w)
% W - Vector of frequencies at which A(w) and B(w) are computed
%
% Description:
%
% The function first prompts the user to select the range of freuencies to
% for idnetification. This is done by clicking with the mouse on a
% plot of the added mass and damping. The low-frequency data point should be clicked first, 
% then the high-frequency point, and finally press return. The selected data is
% then re-plotted.
%
% After selecting the freuquency range, the function allows the elimination
% of data wildpoints. This is done by clicking the points that are to be
% eliminated. Click on all the points to be eliminated and then press
% enter. This process can be re-initialised within the function in case 
% a point in deleted accidentaly.

% You can test this for elements A(6,6) and B(6,6) using the WAMIT tanker data:
%
% > load tanker;
% > EditAB(reshape(vessel.A(6,6,:),60,1),reshape(vessel.B(6,6,:),60,1),vessel.freqs)
%
% Author:    Tristan Perez
% Date:      2009-09-01 
% Revisions: 2021-03-06, minor bug fixes 

% Created by Tristan Perez 
% Date 2009-01-01 
% Revisions:

% Select Frequency range from plot
while 1

    figure(101)
    subplot(211),plot(W,A,'o-','LineWidth',2)
    grid on
    ylabel('A(w)')
    xlabel('Frequency [rad/s]')
    subplot(212),plot(W,B,'o-','LineWidth',2)
    grid on
    ylabel('B(w)')
    xlabel('Frequency [rad/s]')
    disp('---------------------------------------------------------------------')
    disp('Click on the plot twice to define low and high frequency, RETURN when finsihed')
    FR=ginput;
    
    Wmin = FR(1,1);
    Wmax = FR(2,1);
    Dmin=abs(Wmin-W(1));
    Dmax=abs(Wmax-W(1));
    for k=1:length(W)
        Dm = abs(Wmin-W(k));
        DM = abs(Wmax-W(k));
        if Dm<=Dmin
            Dmin=Dm;
            Kmin=k;
        end
        if DM<=Dmax
            Dmax=DM;
            Kmax=k;
        end
    end
    % Crop variables
    if W(Kmin)==0
        Kmin=Kmin+1; %Avoid zero frequency 
    end
    Wm=W(Kmin:Kmax);
    Am=A(Kmin:Kmax);
    Bm=B(Kmin:Kmax);
    
    % Plot added mass and damping to be used for sysid
    figure(101)
    subplot(211),plot(Wm,Am,'o-','LineWidth',2)
    grid on
    ylabel('A(w)')
    xlabel('Frequency [rad/s]')
    subplot(212),plot(Wm,Bm,'o-','LineWidth',2)
    grid on
    ylabel('B(w)')
    xlabel('Frequency [rad/s]')
    
    %% Eliminate wild points
    disp('---------------------------------------------------------------------')
    Usropt=input('Do you need to eliminate wild points? (y = yes or RETURN to continue):','s');
    if Usropt == 'y'
        disp('Click on the wild points and finish with ENTER')
        FR=ginput;
        %
        Wwp = FR(:,1);
        for k=1:length(Wwp)
            for l=1:length(Wm)-1
                if Wm(l) < Wwp(k) && Wwp(k) < Wm(l+1)
                    if abs(Wm(l)-Wwp(k)) < abs(Wwp(k)-Wm(l+1))
                        Wm=[Wm(1:l-1); Wm(l+1:end)];
                        Bm=[Bm(1:l-1); Bm(l+1:end)];
                        Am=[Am(1:l-1); Am(l+1:end)];
                    else
                        Wm=[Wm(1:l); Wm(l+2:end)];
                        Bm=[Bm(1:l); Bm(l+2:end)];
                        Am=[Am(1:l); Am(l+2:end)];
                    end
                end
            end
        end
    else
        break
    end
    % Plot added mass and damping to be used for sysid
    figure(101)
    subplot(211),plot(Wm,Am,'o-','LineWidth',2)
    grid on
    ylabel('A(w)')
    xlabel('Frequency [rad/s]')
    subplot(212),plot(Wm,Bm,'o-','LineWidth',2)
    grid on
    ylabel('B(w)')
    xlabel('Frequency [rad/s]')
    disp('---------------------------------------------------------------------')
    Usropt=input('Start again with range selection an wild point elimination? (y = yes or RETURN to continue):','s');
    if isempty(Usropt)
        break
    end
    
end
    
%% Return Data    
Wmod = Wm;
Amod = Am;
Bmod = Bm;
